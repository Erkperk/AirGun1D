/*
 * airgun.c -- 1D Airgun-Bubble coupled solver
 *
 * SBP-SAT spatial discretization (3rd-order upwind) with RK4 time integration.
 * Real-gas thermodynamics via CoolProp-generated tables + DAK Z-factor.
 * Use --ideal flag for ideal-gas mode (verification against MATLAB code).
 *
 * Usage:
 *   ./airgun                          Real-gas methane (default)
 *   ./airgun --ideal                  Ideal-gas methane (gamma=1.31)
 *   ./airgun --ideal --gamma 1.4 --rgas 287.06  Ideal-gas air (match MATLAB)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* ========== Forward declarations for table-dependent code ========== */
static int g_ideal_gas = 0;

/* ========== Conditionally include tables ========== */
/* Tables are only needed for real-gas mode. In ideal-gas mode the defines
   are still needed to compile, so provide fallback values. */
#if __has_include("tables.h")
#include "tables.h"
#else
/* Fallback: allow compilation without tables (ideal-gas only) */
#define NT 1
#define NP 1
#define T_MIN 150.0
#define T_STEP 2.0
#define P_MIN 10000.0
#define P_STEP 50000.0
static const double U_table[1][1] = {{0.0}};
static const double Cv_table[1][1] = {{0.0}};
#endif

/* ========== Physical constants ========== */
#define M_CH4     16.043       /* g/mol */
#define R_UNIV    8314.46      /* J/(kmol*K) */
#define PI_VAL    3.14159265358979323846

/* Default methane gas properties */
static double g_Rgas   = R_UNIV / M_CH4;   /* 518.28 J/(kg*K) */
static double g_gamma  = 1.31;             /* ratio of heat capacities for CH4 */
static double g_Cv     = 0.0;              /* computed: Rgas / (gamma - 1) */

/* CoolProp U offset: U_table(T,P_low) ≈ Cv*T + U_offset at low pressure.
   Computed at startup so ideal-gas fallback is continuous with tables. */
static double g_U_offset = 0.0;
#define T_TABLE_LO (T_MIN + 5.0 * T_STEP)   /* safe margin inside table */
#define T_TABLE_HI (T_MIN + (NT - 6) * T_STEP)

/* ========== Runtime parameters ========== */
static int    g_nx       = 100;
static double g_pressure = 2000.0 * 6894.76;  /* 2000 psi in Pa */
static double g_length   = 0.5;               /* airgun length [m] */
static double g_area     = 12.5 * 6.4516e-4;  /* port area [m^2] (12.5 in^2) */
static double g_depth    = 10.0;              /* depth [m] */
static double g_endtime  = 2.0;              /* end time [s] */
static double g_cfl      = 5e-4;             /* CFL number */
static double g_cutoff   = 0.30;             /* airgun cutoff time [s] */

/* Water / ambient */
static double g_rho_inf  = 1000.0;  /* water density [kg/m^3] */
static double g_c_inf    = 1482.0;  /* sound speed in water [m/s] */
static double g_pa       = 1e5;     /* atmospheric pressure [Pa] */
static double g_T_inf    = 288.0;   /* ambient temperature [K] */

/* Derived */
static double g_p_inf;     /* ambient pressure at depth */

/* Output */
static int    g_n_output  = 2000;   /* number of output samples */
static double g_r_obs     = 1.0;    /* observation distance [m] */

/* Maximum grid size */
#define MAX_NX 10000

/* ========== Table lookup: bicubic Catmull-Rom interpolation ========== */

static inline double catmull_rom(double fm1, double f0, double f1, double f2, double t)
{
    double t2 = t * t, t3 = t2 * t;
    return 0.5 * ((2.0 * f0) +
                   (-fm1 + f1) * t +
                   (2.0 * fm1 - 5.0 * f0 + 4.0 * f1 - f2) * t2 +
                   (-fm1 + 3.0 * f0 - 3.0 * f1 + f2) * t3);
}

static inline double table_lookup(const double tab[][NP], double T, double P)
{
    if (T < T_MIN) T = T_MIN;
    if (T > T_MIN + (NT - 1) * T_STEP) T = T_MIN + (NT - 1) * T_STEP;
    if (P < P_MIN) P = P_MIN;
    if (P > P_MIN + (NP - 1) * P_STEP) P = P_MIN + (NP - 1) * P_STEP;

    double fi = (T - T_MIN) / T_STEP;
    double fj = (P - P_MIN) / P_STEP;

    int i = (int)fi;
    int j = (int)fj;

    if (i >= NT - 1) i = NT - 2;
    if (j >= NP - 1) j = NP - 2;

    double wt = fi - i;
    double wp = fj - j;

    int im1 = (i > 0) ? i - 1 : 0;
    int ip2 = (i < NT - 2) ? i + 2 : NT - 1;
    int jm1 = (j > 0) ? j - 1 : 0;
    int jp2 = (j < NP - 2) ? j + 2 : NP - 1;

    double r0 = catmull_rom(tab[im1][jm1], tab[im1][j], tab[im1][j+1], tab[im1][jp2], wp);
    double r1 = catmull_rom(tab[i  ][jm1], tab[i  ][j], tab[i  ][j+1], tab[i  ][jp2], wp);
    double r2 = catmull_rom(tab[i+1][jm1], tab[i+1][j], tab[i+1][j+1], tab[i+1][jp2], wp);
    double r3 = catmull_rom(tab[ip2][jm1], tab[ip2][j], tab[ip2][j+1], tab[ip2][jp2], wp);

    return catmull_rom(r0, r1, r2, r3, wt);
}

static inline double lookup_U(double T, double P)
{
    if (g_ideal_gas) return g_Cv * T;
    return table_lookup(U_table, T, P);
}

static inline double lookup_Cv(double T, double P)
{
    if (g_ideal_gas) return g_Cv;
    double cv = table_lookup(Cv_table, T, P);
    if (cv < 500.0) cv = 500.0;
    return cv;
}

/* ========== DAK compression factor (19 coefficients) ========== */

static inline double compression_factor(double T, double P)
{
    if (g_ideal_gas) return 1.0;
    if (T < 150.0) T = 150.0;
    if (P < 100.0) P = 100.0;

    /* Pure methane: gammag = M_CH4 / 28.97 = 0.5538 */
    const double gammag = M_CH4 / 28.97;
    const double Pc = 756.8 - 131.0 * gammag - 3.6 * gammag * gammag;
    const double Tc = 169.2 + 349.5 * gammag - 74.0 * gammag * gammag;

    const double a1  =  0.317842;
    const double a2  =  0.382216;
    const double a3  = -7.76835;
    const double a4  = 14.2905;
    const double a5  =  2e-6;
    const double a6  = -0.004693;
    const double a7  =  0.096254;
    const double a8  =  0.16672;
    const double a9  =  0.96691;
    const double a10 =  0.063069;
    const double a11 = -1.96685;
    const double a12 = 21.0581;
    const double a13 = -27.0246;
    const double a14 = 16.23;
    const double a15 = 207.783;
    const double a16 = -488.161;
    const double a17 = 176.29;
    const double a18 =  1.88453;
    const double a19 =  3.05921;

    double Tpr = (T * 1.8) / Tc;
    double Ppr = (P / 100000.0 * 14.5) / Pc;

    double t = 1.0 / Tpr;

    double A  = a1 * exp(a2 * (1.0 - t) * (1.0 - t)) * Ppr * t;
    double B  = a3 * t + a4 * t * t + a5 * pow(Ppr, 6.0) * pow(t, 6.0);
    double Cc = a9 + a8 * Ppr * t + a7 * Ppr * Ppr * t * t + a6 * Ppr * Ppr * Ppr * t * t * t;
    double Dc = a10 * exp(a11 * (1.0 - t) * (1.0 - t)) * t;
    double Ec = a12 * t + a13 * t * t + a14 * t * t * t;
    double F  = a15 * t + a16 * t * t + a17 * t * t * t;
    double G  = a18 + a19 * t;

    double y = (Dc * Ppr) / (-(A * A * B) / (Cc * Cc * Cc) + (1.0 + A * A) / Cc);

    double Z = (Dc * Ppr * (1.0 + y + y * y - y * y * y))
             / (pow(1.0 - y, 3.0) * (Dc * Ppr + Ec * y * y - F * pow(y, G)));

    return Z;
}

/* ========== Thermodynamic sound speed ========== */

static inline double thermo_soundspeed(double T, double P, double rho_val, double cv)
{
    double gamma_i = 1.0 + g_Rgas / cv;
    if (g_ideal_gas) return sqrt(gamma_i * P / rho_val);
    if (P < 1000.0) return sqrt(gamma_i * P / rho_val);

    double Z = compression_factor(T, P);

    double hT = 0.02;
    double hP = 200.0;
    double dZdT = (compression_factor(T + hT, P) - compression_factor(T - hT, P)) / (2.0 * hT);
    double dZdP = (compression_factor(T, P + hP) - compression_factor(T, P - hP)) / (2.0 * hP);

    double D = 1.0 - P * dZdP / Z;
    if (fabs(D) < 1e-10) D = 1e-10;

    double dpdrho_T = Z * g_Rgas * T / D;
    double dpdT_rho = rho_val * g_Rgas * (T * dZdT + Z) / D;

    double a2 = dpdrho_T + T * dpdT_rho * dpdT_rho / (rho_val * rho_val * cv);

    if (a2 < 0.0) a2 = g_gamma * P / rho_val;
    return sqrt(a2);
}

/* ========== Temperature recovery: Newton iteration on U(T,P) tables ========== */

/* Returns: recovered T, and sets *out_of_table = 1 if ideal-gas fallback was used */
static inline double recover_temperature(double e_int, double T_guess, double P_guess,
                                          int *out_of_table)
{
    *out_of_table = 0;
    if (g_ideal_gas) return e_int / g_Cv;

    double T = T_guess;
    double P = (P_guess < P_MIN) ? P_MIN : P_guess;

    for (int iter = 0; iter < 5; iter++) {
        double cv = lookup_Cv(T, P);
        T += (e_int - lookup_U(T, P)) / cv;
        if (T < T_MIN) T = T_MIN;
        if (T > T_MIN + (NT - 1) * T_STEP) T = T_MIN + (NT - 1) * T_STEP;
    }

    /* If T lands outside safe table interior, fall back to ideal gas */
    if (T < T_TABLE_LO || T > T_TABLE_HI || T != T) {
        T = (e_int - g_U_offset) / g_Cv;
        if (T < 50.0) T = 50.0;
        *out_of_table = 1;
    }

    return T;
}

/* ========== Pressure from state ========== */

static inline double compute_pressure(double rho, double T)
{
    if (g_ideal_gas) return rho * g_Rgas * T;
    double Z = compression_factor(T, fmax(rho * g_Rgas * T, P_MIN));
    /* Iterate once to get consistent P */
    double P = rho * Z * g_Rgas * T;
    Z = compression_factor(T, fmax(P, P_MIN));
    P = rho * Z * g_Rgas * T;
    return P;
}

/* ========== SBP Operators (3rd-order upwind) ========== */
/*
 * From sbplib/+sbp/+implementations/d1_upwind_3.m
 *
 * H norm diagonal: boundary [3/8, 7/6, 23/24], interior 1, times h
 * Qp boundary block (3x3):
 *   [-1/24,  17/24,  -1/6 ]
 *   [-13/24, -1/4,   23/24]
 *   [ 1/12,  -11/24, -11/24]
 * Qp interior stencil (offsets -1,0,+1,+2): [-1/3, -1/2, 1, -1/6]
 * Qm = -Qp^T
 *
 * Dp = H^{-1} * (Qp - 1/2*e1*e1' + 1/2*em*em')
 * Dm = H^{-1} * (Qm - 1/2*e1*e1' + 1/2*em*em')
 *
 * D1 = (Dp + Dm) / 2  (central)
 * Ddisp = (Dp - Dm) / 2  (dissipation)
 */

/* H norm diagonal (stored as 1/h factors; actual H = Hv[i] * h) */
static double g_Hv[MAX_NX];
static double g_HIv[MAX_NX];  /* H inverse: 1/(Hv[i]*h) */
static double g_h;             /* grid spacing */
static int    g_m;             /* number of grid points */

static void init_sbp(int m, double h)
{
    g_m = m;
    g_h = h;

    /* H norm diagonal */
    for (int i = 0; i < m; i++) g_Hv[i] = 1.0;

    g_Hv[0]   = 3.0 / 8.0;
    g_Hv[1]   = 7.0 / 6.0;
    g_Hv[2]   = 23.0 / 24.0;
    g_Hv[m-3] = 23.0 / 24.0;
    g_Hv[m-2] = 7.0 / 6.0;
    g_Hv[m-1] = 3.0 / 8.0;

    for (int i = 0; i < m; i++) g_HIv[i] = 1.0 / (g_Hv[i] * h);
}

/*
 * Apply Qp to a vector f of length m, result in out.
 *
 * Built as: spdiags fills interior stencil [-1/3,-1/2,1,-1/6] at offsets
 * [-1,0,+1,+2] for ALL rows, then boundary 3x3 blocks overwrite columns
 * 0..2 of rows 0..2 and columns m-3..m-1 of rows m-3..m-1.
 * This leaves interior stencil entries BEYOND the 3x3 blocks intact.
 *
 * Full row stencils (verified: all row sums = 0 except first/last = ±1/2):
 *   Row 0:   [-1/24, 17/24, -1/6]                    (cols 0..2)
 *   Row 1:   [-13/24, -1/4, 23/24, -1/6]             (cols 0..3)
 *   Row 2:   [1/12, -11/24, -11/24, 1, -1/6]         (cols 0..4)
 *   Row i:   [-1/3, -1/2, 1, -1/6]                   (cols i-1..i+2)
 *   Row m-3: [-1/3, -11/24, 23/24, -1/6]             (cols m-4..m-1)
 *   Row m-2: [-11/24, -1/4, 17/24]                   (cols m-3..m-1)
 *   Row m-1: [1/12, -13/24, -1/24]                   (cols m-3..m-1)
 */
static void apply_Qp(const double *f, double *out, int m)
{
    /* Top boundary rows */
    out[0] = (-1.0/24.0)*f[0] + (17.0/24.0)*f[1] + (-1.0/6.0)*f[2];

    out[1] = (-13.0/24.0)*f[0] + (-1.0/4.0)*f[1] + (23.0/24.0)*f[2];
    if (3 < m) out[1] += (-1.0/6.0)*f[3];

    out[2] = (1.0/12.0)*f[0] + (-11.0/24.0)*f[1] + (-11.0/24.0)*f[2];
    if (3 < m) out[2] += (1.0)*f[3];
    if (4 < m) out[2] += (-1.0/6.0)*f[4];

    /* Interior rows */
    for (int i = 3; i <= m - 4; i++) {
        out[i] = (-1.0/3.0)*f[i-1] + (-1.0/2.0)*f[i] + (1.0)*f[i+1] + (-1.0/6.0)*f[i+2];
    }

    /* Bottom boundary rows */
    out[m-3] = (-11.0/24.0)*f[m-3] + (23.0/24.0)*f[m-2] + (-1.0/6.0)*f[m-1];
    if (m-4 >= 0) out[m-3] += (-1.0/3.0)*f[m-4];

    out[m-2] = (-11.0/24.0)*f[m-3] + (-1.0/4.0)*f[m-2] + (17.0/24.0)*f[m-1];

    out[m-1] = (1.0/12.0)*f[m-3] + (-13.0/24.0)*f[m-2] + (-1.0/24.0)*f[m-1];
}

/*
 * Apply Qm = -Qp^T to a vector f of length m.
 * Interior stencil at offsets [-2,-1,0,+1]: [1/6, -1, 1/2, 1/3]
 *
 * Full row stencils:
 *   Row 0:   [1/24, 13/24, -1/12]                       (cols 0..2)
 *   Row 1:   [-17/24, 1/4, 11/24]                       (cols 0..2)
 *   Row 2:   [1/6, -23/24, 11/24, 1/3]                  (cols 0..3)
 *   Row j:   [1/6, -1, 1/2, 1/3]                        (cols j-2..j+1)
 *   Row m-3: [1/6, -1, 11/24, 11/24, -1/12]             (cols m-5..m-1)
 *   Row m-2: [1/6, -23/24, 1/4, 13/24]                  (cols m-4..m-1)
 *   Row m-1: [1/6, -17/24, 1/24]                        (cols m-3..m-1)
 */
static void apply_Qm(const double *f, double *out, int m)
{
    /* Top boundary rows */
    out[0] = (1.0/24.0)*f[0] + (13.0/24.0)*f[1] + (-1.0/12.0)*f[2];

    out[1] = (-17.0/24.0)*f[0] + (1.0/4.0)*f[1] + (11.0/24.0)*f[2];

    out[2] = (1.0/6.0)*f[0] + (-23.0/24.0)*f[1] + (11.0/24.0)*f[2];
    if (3 < m) out[2] += (1.0/3.0)*f[3];

    /* Interior rows */
    for (int j = 3; j <= m - 4; j++) {
        out[j] = (1.0/6.0)*f[j-2] + (-1.0)*f[j-1] + (1.0/2.0)*f[j] + (1.0/3.0)*f[j+1];
    }

    /* Bottom boundary rows */
    out[m-3] = (11.0/24.0)*f[m-3] + (11.0/24.0)*f[m-2] + (-1.0/12.0)*f[m-1];
    if (m-4 >= 0) out[m-3] += (-1.0)*f[m-4];
    if (m-5 >= 0) out[m-3] += (1.0/6.0)*f[m-5];

    out[m-2] = (-23.0/24.0)*f[m-3] + (1.0/4.0)*f[m-2] + (13.0/24.0)*f[m-1];
    if (m-4 >= 0) out[m-2] += (1.0/6.0)*f[m-4];

    out[m-1] = (1.0/6.0)*f[m-3] + (-17.0/24.0)*f[m-2] + (1.0/24.0)*f[m-1];
}

/*
 * Apply Dp = H^{-1} * (Qp + SAT) to interleaved data q[3*m].
 * SAT boundary terms: -1/2 at (0,0), +1/2 at (m-1,m-1).
 * Operates via Kronecker D ⊗ I_3: each component independently.
 */
static void apply_Dp(const double *q, double *out, int m)
{
    double f[MAX_NX], r[MAX_NX];
    for (int c = 0; c < 3; c++) {
        /* Extract component c */
        for (int i = 0; i < m; i++) f[i] = q[3*i + c];

        /* Apply Qp */
        apply_Qp(f, r, m);

        /* Add SAT: -1/2*e1*e1' + 1/2*em*em' */
        r[0]   -= 0.5 * f[0];
        r[m-1] += 0.5 * f[m-1];

        /* Apply H^{-1} */
        for (int i = 0; i < m; i++) out[3*i + c] = g_HIv[i] * r[i];
    }
}

static void apply_Dm(const double *q, double *out, int m)
{
    double f[MAX_NX], r[MAX_NX];
    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < m; i++) f[i] = q[3*i + c];

        apply_Qm(f, r, m);

        r[0]   -= 0.5 * f[0];
        r[m-1] += 0.5 * f[m-1];

        for (int i = 0; i < m; i++) out[3*i + c] = g_HIv[i] * r[i];
    }
}

/*
 * D1 = (Dp + Dm) / 2   -- central derivative
 * Ddisp = (Dp - Dm) / 2  -- dissipation operator
 *
 * Euler flux with upwind dissipation:
 *   dq = -D1 * F(q) + Ddisp * (lambda_max * q)
 */

/* Working arrays for the Euler flux computation */
static double g_flux[3 * MAX_NX];    /* F(q) interleaved */
static double g_lamb_q[3 * MAX_NX];  /* lambda_max * q interleaved */
/* (removed unused work arrays) */

/* Per-grid-point thermodynamic state (cached for BC computation) */
static double g_T[MAX_NX];
static double g_p[MAX_NX];
static double g_c[MAX_NX];

static void compute_euler_rhs(const double *q, double *dq, int m)
{
    /* Recover thermodynamic state and build flux/dissipation vectors */
    for (int i = 0; i < m; i++) {
        double rho  = q[3*i + 0];
        double rhou = q[3*i + 1];
        double e    = q[3*i + 2];

        if (rho < 1e-10) rho = 1e-10;

        double u = rhou / rho;
        double e_int = e / rho - 0.5 * u * u;
        if (e_int < 100.0) e_int = 100.0;

        /* Recover temperature */
        double T_guess = g_T[i] > 0.0 ? g_T[i] : g_T_inf;
        double P_guess = g_p[i] > 0.0 ? g_p[i] : g_pressure;
        int oot = 0;
        double T = recover_temperature(e_int, T_guess, P_guess, &oot);
        double P, cv, c;
        if (oot) {
            /* Out-of-table fallback: ideal-gas EOS */
            P = rho * g_Rgas * T;
            cv = g_Cv;
            c = sqrt(g_gamma * P / rho);
        } else {
            P = compute_pressure(rho, T);
            cv = lookup_Cv(T, fmax(P, P_MIN));
            c = thermo_soundspeed(T, P, rho, cv);
        }

        g_T[i] = T;
        g_p[i] = P;
        g_c[i] = c;

        /* Euler flux: F = [rho*u, rho*u^2 + p, (e+p)*u] */
        g_flux[3*i + 0] = rhou;
        g_flux[3*i + 1] = rhou * u + P;
        g_flux[3*i + 2] = (e + P) * u;

        /* Dissipation: lambda_max * q */
        double lambda_max = c + fabs(u);
        g_lamb_q[3*i + 0] = lambda_max * rho;
        g_lamb_q[3*i + 1] = lambda_max * rhou;
        g_lamb_q[3*i + 2] = lambda_max * e;
    }

    /* Apply SBP operators: dq = -D1*F + Ddisp*(lambda*q)
     *   D1 = (Dp + Dm)/2, Ddisp = (Dp - Dm)/2
     *   => dq = -(Dp+Dm)/2 * F + (Dp-Dm)/2 * (lambda*q)
     *        = (-Dp*F + Dp*(lambda*q) - Dm*F - Dm*(lambda*q)) / 2   NO
     *   Actually: dq = -0.5*(Dp+Dm)*F + 0.5*(Dp-Dm)*(lambda*q)
     *              = 0.5*(-Dp*F - Dm*F + Dp*(lambda*q) - Dm*(lambda*q))
     *              = 0.5*(Dp*(-F + lambda*q) + Dm*(-F - lambda*q))
     * But simpler to compute Dp*F, Dm*F, Dp*(lambda*q), Dm*(lambda*q) and combine.
     * Actually even simpler: just compute Dp and Dm on F and lamb_q separately.
     */

    /* We need: -D1*F + Ddisp*lamb_q = -0.5*(Dp+Dm)*F + 0.5*(Dp-Dm)*lamb_q
     *        = 0.5*Dp*(-F+lamb_q) + 0.5*Dm*(-F-lamb_q) ... no, let me just be explicit. */

    /* Compute Dp*F and Dm*F */
    double DpF[3 * MAX_NX], DmF[3 * MAX_NX];
    apply_Dp(g_flux, DpF, m);
    apply_Dm(g_flux, DmF, m);

    /* Compute Dp*(lambda*q) and Dm*(lambda*q) */
    double DpLQ[3 * MAX_NX], DmLQ[3 * MAX_NX];
    apply_Dp(g_lamb_q, DpLQ, m);
    apply_Dm(g_lamb_q, DmLQ, m);

    /* dq = -0.5*(DpF + DmF) + 0.5*(DpLQ - DmLQ) */
    for (int k = 0; k < 3 * m; k++) {
        dq[k] = -0.5 * (DpF[k] + DmF[k]) + 0.5 * (DpLQ[k] - DmLQ[k]);
    }
}

/* ========== Characteristic SAT boundary conditions ========== */

/* Effective gamma for real gas: gamma_eff = p / (e - 0.5*rho*u^2) + 1 */
static inline double effective_gamma(double rho, double rhou, double e, double p)
{
    double e_int_total = e - 0.5 * rhou * rhou / rho;
    if (e_int_total < 1e-6) e_int_total = 1e-6;
    return p / e_int_total + 1.0;
}

/*
 * Left wall BC (x = -L, u = 0): Reflect incoming characteristic.
 *
 * From Euler1d.m boundary_condition_wall:
 *   At left boundary (s = -1):
 *     p_in = 2 (eigenvalue u+c, incoming from left)
 *     p_ot = 3 (eigenvalue u-c, outgoing)
 *     p_zero = 1 (entropy wave, eigenvalue u)
 *   T matrix columns: eigenvectors of flux Jacobian
 *   R = -(u-c)/(u+c)  -- reflection coefficient
 *   tau1 = -2*c (penalty on incoming characteristic)
 *   SAT = 1/2 * H^{-1} * tau * (w_in - R*w_out)
 *
 *   where tau = -s * e_S * T * tauHat(pt) = e_S * T * [-2c; 0; 0] permuted
 */
static void apply_wall_bc_left(const double *q, double *sat, int m)
{
    (void)m;
    double rho  = q[0];
    double rhou = q[1];
    double e    = q[2];

    if (rho < 1e-10) rho = 1e-10;
    double u = rhou / rho;
    double gamma_e = effective_gamma(rho, rhou, e, g_p[0]);
    double c = g_c[0];

    double sqrt2gm = sqrt(2.0 * (gamma_e - 1.0));

    /* T matrix (columns are eigenvectors for eigenvalues u, u+c, u-c) */
    double T11 = sqrt2gm * rho;
    double T12 = rho;
    double T13 = rho;
    double T21 = sqrt2gm * rho * u;
    double T22 = rho * (u + c);
    double T23 = rho * (u - c);
    double T31 = sqrt2gm * rho * u * u / 2.0;
    double T32 = e + (gamma_e - 1.0) * (e - rho * u * u / 2.0) + rho * u * c;
    double T33 = e + (gamma_e - 1.0) * (e - rho * u * u / 2.0) - rho * u * c;

    /* T^{-1} via Cramer's rule (3x3) */
    double det = T11*(T22*T33 - T23*T32) - T12*(T21*T33 - T23*T31) + T13*(T21*T32 - T22*T31);
    if (fabs(det) < 1e-30) det = 1e-30;
    double idet = 1.0 / det;

    /* Only rows 2,3 of T^{-1} needed (w2, w3 characteristics) */
    double Ti21 = (T23*T31 - T21*T33) * idet;
    double Ti22 = (T11*T33 - T13*T31) * idet;
    double Ti23 = (T13*T21 - T11*T23) * idet;
    double Ti31 = (T21*T32 - T22*T31) * idet;
    double Ti32 = (T12*T31 - T11*T32) * idet;
    double Ti33 = (T11*T22 - T12*T21) * idet;

    /* Characteristic variables w = T^{-1} * q */
    double w2 = Ti21*rho + Ti22*rhou + Ti23*e;
    double w3 = Ti31*rho + Ti32*rhou + Ti33*e;

    /* At left boundary (s=-1): incoming = w2 (u+c > 0), outgoing = w3 (u-c < 0)
     * Permutation: p = [p_in=2, p_zero=1, p_ot=3]
     * R = -(u-c)/(u+c) */
    double R_coeff = -(u - c) / (u + c);
    double w_in = w2;
    double w_ot = w3;

    /* Penalty tauHat = [-2c, 0, 0] permuted back:
     * tauHat[p_in] = -2c, tauHat[p_zero] = 0, tauHat[p_ot] = 0
     * After inverse permutation (pt): tauHat_orig[1] = 0, tauHat_orig[2] = -2c, tauHat_orig[3] = 0
     * tau_vec = -s * T * tauHat = 1 * T * [0, -2c, 0]^T = -2c * T_col2
     *
     * Actually s = -1 for left boundary.
     * tau = -s * e_S * T * tauHat(pt)
     * tauHat = [tau1=(-2c); tau2=[0;0]] in permuted space
     * tauHat(pt) maps back: position p_in=2 gets -2c, others get 0
     * So tauHat_orig = [0; -2*c; 0]
     * tau = -(-1) * T * [0; -2*c; 0] = T * [0; -2*c; 0] = -2*c * [T12; T22; T32]
     */
    double tau0 = -2.0 * c * T12;
    double tau1 = -2.0 * c * T22;
    double tau2 = -2.0 * c * T32;

    /* SAT = 1/2 * H^{-1} * tau * (w_in - R * w_ot) */
    double strength = 0.5 * (w_in - R_coeff * w_ot);
    double HI0 = g_HIv[0];

    sat[0] = HI0 * tau0 * strength;
    sat[1] = HI0 * tau1 * strength;
    sat[2] = HI0 * tau2 * strength;
}

/*
 * Right outflow BC (x = 0, set bubble pressure).
 *
 * From Euler1d.m boundary_condition_outflow -> boundary_condition_L:
 *   At right boundary (s = 1):
 *     p_in = 3 (eigenvalue u-c, incoming from right when u < c)
 *   L = (gamma-1) * [0, -u/2, 1]  (extracts pressure)
 *   tau = e_S * T * tauHat(pt)  where tauHat = [-2|Lambda_in|; 0; 0]
 *   SAT = 1/2 * H^{-1} * tau * inv(L*T_in) * (L*q - p_bubble)
 */

/* Flow state at right boundary */
#define SUBSONIC_INFLOW    1
#define SUBSONIC_OUTFLOW  (-1)
#define SUPERSONIC_OUTFLOW (-2)

static int flow_state_right(const double *q, int m)
{
    double rho  = q[3*(m-1) + 0];
    double rhou = q[3*(m-1) + 1];

    if (rho < 1e-10) rho = 1e-10;
    double u = rhou / rho;
    double c = g_c[m-1];

    if (u < -c) return 2;   /* supersonic inflow */
    if (u < 0)  return SUBSONIC_INFLOW;
    if (u < c)  return SUBSONIC_OUTFLOW;
    return SUPERSONIC_OUTFLOW;
}

static void apply_outflow_bc_right(const double *q, double *sat, double p_bubble, int m)
{
    int idx = 3 * (m - 1);
    double rho  = q[idx + 0];
    double rhou = q[idx + 1];
    double e    = q[idx + 2];

    if (rho < 1e-10) rho = 1e-10;
    double u = rhou / rho;
    double gamma_e = effective_gamma(rho, rhou, e, g_p[m-1]);
    double c = g_c[m-1];

    /* T matrix column 3 (eigenvector for u-c, the incoming characteristic) */
    double T13 = rho;
    double T23 = rho * (u - c);
    double T33 = e + (gamma_e - 1.0) * (e - rho * u * u / 2.0) - rho * u * c;

    /* At right boundary (s=1): incoming char = 3 (u-c, which is < 0 for subsonic outflow)
     * p_in = 3, p_ot = [1, 2]
     * Permutation: p = [3, 1, 2], inverse pt: pt[3]=1, pt[1]=2, pt[2]=3
     * tauHat = [-2*|u-c|; 0; 0] in permuted space
     * tauHat(pt): position 3 gets -2|u-c|, others 0
     * => tauHat_orig = [0; 0; -2*|u-c|]
     * tau = e_S * T * tauHat_orig = -2*|u-c| * T_col3
     */
    double abs_lambda_in = fabs(u - c);
    double tau0 = -2.0 * abs_lambda_in * T13;
    double tau1 = -2.0 * abs_lambda_in * T23;
    double tau2 = -2.0 * abs_lambda_in * T33;

    /* L = (gamma_e - 1) * [0, -u/2, 1] */
    double gm1 = gamma_e - 1.0;
    double L0 = 0.0;
    double L1 = gm1 * (-u / 2.0);
    double L2 = gm1;

    /* L*T_in = L * T_col3 */
    double LTin = L0 * T13 + L1 * T23 + L2 * T33;
    if (fabs(LTin) < 1e-30) LTin = 1e-30;

    /* L*q_boundary */
    double Lq = L0 * rho + L1 * rhou + L2 * e;

    /* SAT = 1/2 * H^{-1} * tau * (1/LTin) * (Lq - p_bubble) */
    double strength = 0.5 * (Lq - p_bubble) / LTin;
    double HI_end = g_HIv[m-1];

    sat[idx + 0] = HI_end * tau0 * strength;
    sat[idx + 1] = HI_end * tau1 * strength;
    sat[idx + 2] = HI_end * tau2 * strength;
}

/* ========== Bubble dynamics ========== */
/*
 * State: [R, Rdot, m, E]
 * From bubbleRHS.m and bubblePressure.m
 */

static inline double bubble_pressure(const double *bubble)
{
    double R = bubble[0];
    double E = bubble[3];
    double V = (4.0 / 3.0) * PI_VAL * R * R * R;
    return E * (g_gamma - 1.0) / V;
}

static void bubble_rhs(const double *bubble, double rho_a, double v_a,
                        double e_a, double p_a, double A,
                        double *dbubble)
{
    double R    = bubble[0];
    double Rdot = bubble[1];
    double m_b  = bubble[2];
    double E    = bubble[3];

    double p_b   = bubble_pressure(bubble);
    double V     = (4.0 / 3.0) * PI_VAL * R * R * R;
    double Vdot  = 4.0 * PI_VAL * R * R * Rdot;

    /* Heat transfer */
    double kappa = 4000.0;
    double M_ht  = 10.0;
    double T_water = 273.0;
    double cv_bubble = 718.0;  /* bubble gas Cv (simple model) */
    double Tb = E / (cv_bubble * m_b);
    double dQdt = 4.0 * PI_VAL * R * R * M_ht * kappa * (Tb - T_water);

    /* Rates */
    double dR = Rdot;
    double dE = A * (e_a + p_a) * v_a - p_b * Vdot - dQdt;

    double dpdt = (g_gamma - 1.0) * (dE * V - Vdot * E) / (V * V);

    /* Keller-Miksis with damping */
    double alpha = 0.8;
    double dRdot = (1.0 / R) * ((p_b - g_p_inf) / g_rho_inf
                   + R / (g_rho_inf * g_c_inf) * dpdt
                   - 1.5 * Rdot * Rdot
                   - alpha * Rdot);

    double dm = A * rho_a * v_a;

    dbubble[0] = dR;
    dbubble[1] = dRdot;
    dbubble[2] = dm;
    dbubble[3] = dE;
}

/* ========== Coupled RHS ========== */
/*
 * State vector: [q(3*m), bubble(4)]
 * Total DOF: 3*m + 4
 */

static void coupled_rhs(const double *state, double t, double *dstate, int m)
{
    int nq = 3 * m;
    const double *q = state;
    const double *bubble = state + nq;

    double *dq = dstate;
    double *dbubble = dstate + nq;

    /* Flow state at right boundary */
    int fs = flow_state_right(q, m);

    /* After cutoff or subsonic inflow: chamber closed */
    if (t >= g_cutoff || fs == SUBSONIC_INFLOW) {
        for (int k = 0; k < nq; k++) dq[k] = 0.0;
        /* Bubble still evolves but with zero mass/energy flux */
        bubble_rhs(bubble, 0.0, 0.0, 0.0, 0.0, 0.0, dbubble);
        return;
    }

    /* Bubble pressure for outflow BC */
    double p_b = bubble_pressure(bubble);

    /* Euler flux + dissipation */
    compute_euler_rhs(q, dq, m);

    /* Wall BC at left */
    double sat_l[3] = {0.0, 0.0, 0.0};
    apply_wall_bc_left(q, sat_l, m);
    dq[0] += sat_l[0];
    dq[1] += sat_l[1];
    dq[2] += sat_l[2];

    /* Outflow BC at right (subsonic outflow) */
    if (fs == SUBSONIC_OUTFLOW) {
        double sat_r[3 * MAX_NX];
        memset(sat_r, 0, nq * sizeof(double));
        apply_outflow_bc_right(q, sat_r, p_b, m);
        for (int k = 0; k < nq; k++) dq[k] += sat_r[k];
    }
    /* Supersonic outflow: no BC needed */

    /* Extract state at right boundary for bubble coupling */
    int idx_r = 3 * (m - 1);
    double rho_a = q[idx_r + 0];
    double rhou_a = q[idx_r + 1];
    double e_a = q[idx_r + 2];
    if (rho_a < 1e-10) rho_a = 1e-10;
    double v_a = rhou_a / rho_a;
    double p_a = g_p[m-1];

    bubble_rhs(bubble, rho_a, v_a, e_a, p_a, g_area, dbubble);
}

/* ========== RK4 time integration ========== */

static double *g_state;    /* current state */
static double *g_k1, *g_k2, *g_k3, *g_k4;  /* RK stages */
static double *g_tmp;      /* temporary state */

static void rk4_step(double t, double dt, int ndof)
{
    int m = g_m;

    /* k1 = f(t, y) */
    coupled_rhs(g_state, t, g_k1, m);

    /* k2 = f(t + dt/2, y + dt/2 * k1) */
    for (int i = 0; i < ndof; i++) g_tmp[i] = g_state[i] + 0.5 * dt * g_k1[i];
    coupled_rhs(g_tmp, t + 0.5 * dt, g_k2, m);

    /* k3 = f(t + dt/2, y + dt/2 * k2) */
    for (int i = 0; i < ndof; i++) g_tmp[i] = g_state[i] + 0.5 * dt * g_k2[i];
    coupled_rhs(g_tmp, t + 0.5 * dt, g_k3, m);

    /* k4 = f(t + dt, y + dt * k3) */
    for (int i = 0; i < ndof; i++) g_tmp[i] = g_state[i] + dt * g_k3[i];
    coupled_rhs(g_tmp, t + dt, g_k4, m);

    /* y_{n+1} = y_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4) */
    for (int i = 0; i < ndof; i++) {
        g_state[i] += (dt / 6.0) * (g_k1[i] + 2.0 * g_k2[i] + 2.0 * g_k3[i] + g_k4[i]);
    }
}

/* ========== Output ========== */

static void print_header(void)
{
    printf("# t  R  Rdot  p_bubble  p_acoustic\n");
}

static void print_state(double t, const double *state, const double *dstate, int m)
{
    int nq = 3 * m;
    const double *bubble = state + nq;
    const double *dbubble = dstate + nq;

    double R    = bubble[0];
    double Rdot = bubble[1];
    double p_b  = bubble_pressure(bubble);

    /* Acoustic pressure: p_ac = rho_inf * (2*R*Rdot^2 + R^2*Rddot) / r */
    double Rddot = dbubble[1];
    double p_ac = g_rho_inf * (2.0 * R * Rdot * Rdot + R * R * Rddot) / g_r_obs;

    printf("%.8e  %.8e  %.8e  %.8e  %.8e\n", t, R, Rdot, p_b, p_ac);
}

/* ========== CLI parsing ========== */

static void print_usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Options:\n"
        "  --nx <int>        Number of grid points per meter  [100]\n"
        "  --pressure <Pa>   Chamber pressure                 [2000 psi]\n"
        "  --length <m>      Airgun chamber length             [0.5]\n"
        "  --area <m^2>      Port area                         [12.5 in^2]\n"
        "  --depth <m>       Water depth                       [10]\n"
        "  --endtime <s>     Simulation end time               [2.0]\n"
        "  --cfl <val>       CFL number                        [5e-4]\n"
        "  --cutoff <s>      Airgun cutoff time                [0.30]\n"
        "  --ideal           Use ideal-gas EOS\n"
        "  --gamma <val>     Ratio of heat capacities          [1.31]\n"
        "  --rgas <val>      Specific gas constant [J/(kg*K)]  [518.28]\n"
        "  --robs <m>        Observation distance               [1.0]\n"
        "  --nout <int>      Number of output samples           [2000]\n"
        "\n", prog);
}

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--nx") == 0 && i + 1 < argc) {
            g_nx = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--pressure") == 0 && i + 1 < argc) {
            g_pressure = atof(argv[++i]);
        } else if (strcmp(argv[i], "--length") == 0 && i + 1 < argc) {
            g_length = atof(argv[++i]);
        } else if (strcmp(argv[i], "--area") == 0 && i + 1 < argc) {
            g_area = atof(argv[++i]);
        } else if (strcmp(argv[i], "--depth") == 0 && i + 1 < argc) {
            g_depth = atof(argv[++i]);
        } else if (strcmp(argv[i], "--endtime") == 0 && i + 1 < argc) {
            g_endtime = atof(argv[++i]);
        } else if (strcmp(argv[i], "--cfl") == 0 && i + 1 < argc) {
            g_cfl = atof(argv[++i]);
        } else if (strcmp(argv[i], "--cutoff") == 0 && i + 1 < argc) {
            g_cutoff = atof(argv[++i]);
        } else if (strcmp(argv[i], "--ideal") == 0) {
            g_ideal_gas = 1;
        } else if (strcmp(argv[i], "--gamma") == 0 && i + 1 < argc) {
            g_gamma = atof(argv[++i]);
        } else if (strcmp(argv[i], "--rgas") == 0 && i + 1 < argc) {
            g_Rgas = atof(argv[++i]);
        } else if (strcmp(argv[i], "--robs") == 0 && i + 1 < argc) {
            g_r_obs = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nout") == 0 && i + 1 < argc) {
            g_n_output = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            exit(0);
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            exit(1);
        }
    }
}

/* ========== Main ========== */

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    /* Derived constants */
    g_Cv = g_Rgas / (g_gamma - 1.0);
    g_p_inf = g_pa + g_rho_inf * 9.8 * g_depth;

    /* Grid: m points on [-L, 0] */
    int m = (int)(g_length * g_nx) + 1;
    if (m < 6) m = 6;
    if (m > MAX_NX) {
        fprintf(stderr, "Error: m=%d exceeds MAX_NX=%d\n", m, MAX_NX);
        return 1;
    }

    double h = g_length / (m - 1);

    fprintf(stderr, "Airgun 1D solver\n");
    fprintf(stderr, "  Mode: %s\n", g_ideal_gas ? "ideal gas" : "real gas (tables)");
    fprintf(stderr, "  Gas: gamma=%.4f, Rgas=%.2f, Cv=%.2f\n", g_gamma, g_Rgas, g_Cv);
    fprintf(stderr, "  Chamber: P=%.1f bar, L=%.3f m, A=%.6f m^2\n",
            g_pressure / 1e5, g_length, g_area);
    fprintf(stderr, "  Depth: %.1f m, p_inf=%.1f bar\n", g_depth, g_p_inf / 1e5);
    fprintf(stderr, "  Grid: m=%d, h=%.6f m\n", m, h);
    fprintf(stderr, "  CFL: %.2e, cutoff: %.2f s, endtime: %.2f s\n",
            g_cfl, g_cutoff, g_endtime);

    if (!g_ideal_gas) {
        /* Verify tables are loaded (check NT > 1) */
        if (NT <= 1) {
            fprintf(stderr, "Error: tables.h not available. Use --ideal or generate tables.\n");
            return 1;
        }
        /* Compute U offset so ideal-gas fallback is continuous with tables.
         * At low pressure, U(T,P_low) ≈ Cv*T + U_offset. */
        double T_ref_lo = 300.0;
        g_U_offset = lookup_U(T_ref_lo, P_MIN) - g_Cv * T_ref_lo;
        fprintf(stderr, "  U_offset: %.1f J/kg (ideal-gas fallback)\n", g_U_offset);
    }

    /* Initialize SBP operators */
    init_sbp(m, h);

    /* Allocate state: 3*m (Euler) + 4 (bubble) */
    int ndof = 3 * m + 4;
    g_state = (double *)calloc(ndof, sizeof(double));
    g_k1    = (double *)calloc(ndof, sizeof(double));
    g_k2    = (double *)calloc(ndof, sizeof(double));
    g_k3    = (double *)calloc(ndof, sizeof(double));
    g_k4    = (double *)calloc(ndof, sizeof(double));
    g_tmp   = (double *)calloc(ndof, sizeof(double));

    if (!g_state || !g_k1 || !g_k2 || !g_k3 || !g_k4 || !g_tmp) {
        fprintf(stderr, "Error: memory allocation failed\n");
        return 1;
    }

    /* ---- Initial conditions (from configAirgun.m) ---- */

    /* Airgun chamber: uniform p0, T_inf, u=0 */
    double T0 = g_T_inf;
    double rho0, e0;

    if (g_ideal_gas) {
        rho0 = g_pressure / (g_Rgas * T0);
        e0 = g_Cv * rho0 * T0;
    } else {
        double Z0 = compression_factor(T0, g_pressure);
        rho0 = g_pressure / (Z0 * g_Rgas * T0);
        double u_int0 = lookup_U(T0, g_pressure);
        e0 = rho0 * u_int0;
    }

    for (int i = 0; i < m; i++) {
        g_state[3*i + 0] = rho0;
        g_state[3*i + 1] = 0.0;      /* rho*u = 0 */
        g_state[3*i + 2] = e0;
    }

    /* Initialize thermodynamic cache */
    for (int i = 0; i < m; i++) {
        g_T[i] = T0;
        g_p[i] = g_pressure;
        double cv = g_ideal_gas ? g_Cv : lookup_Cv(T0, g_pressure);
        g_c[i] = g_ideal_gas ? sqrt(g_gamma * g_pressure / rho0)
                             : thermo_soundspeed(T0, g_pressure, rho0, cv);
    }

    fprintf(stderr, "  Initial: rho=%.4f, e=%.1f, c=%.1f m/s\n", rho0, e0, g_c[0]);

    /* Bubble: volume = airgun volume, at ambient conditions */
    double V_gun = g_length * g_area;
    double R_bub = pow(3.0 / (4.0 * PI_VAL) * V_gun, 1.0 / 3.0);
    double m_bub = g_p_inf * V_gun / (g_Rgas * T0);
    double E_bub = g_Cv * m_bub * T0;

    int nq = 3 * m;
    g_state[nq + 0] = R_bub;
    g_state[nq + 1] = 0.0;
    g_state[nq + 2] = m_bub;
    g_state[nq + 3] = E_bub;

    fprintf(stderr, "  Bubble: R=%.6f m, m=%.6f kg, E=%.1f J, p=%.1f bar\n",
            R_bub, m_bub, E_bub, bubble_pressure(g_state + nq) / 1e5);

    /* ---- Time integration ---- */
    double dt = g_cfl * h;
    double t = 0.0;
    int n_steps = (int)(g_endtime / dt) + 1;
    int output_interval = n_steps / g_n_output;
    if (output_interval < 1) output_interval = 1;

    fprintf(stderr, "  dt=%.2e, steps=%d, output every %d steps\n",
            dt, n_steps, output_interval);

    print_header();

    /* Print initial state */
    coupled_rhs(g_state, t, g_k1, m);
    print_state(t, g_state, g_k1, m);

    clock_t t_start = clock();

    for (int step = 1; step <= n_steps; step++) {
        rk4_step(t, dt, ndof);
        t += dt;

        if (step % output_interval == 0 || step == n_steps) {
            coupled_rhs(g_state, t, g_k1, m);
            print_state(t, g_state, g_k1, m);
        }

        /* Progress */
        if (step % (n_steps / 10 + 1) == 0) {
            double elapsed = (double)(clock() - t_start) / CLOCKS_PER_SEC;
            fprintf(stderr, "  step %d/%d (t=%.4f), elapsed %.1f s\n",
                    step, n_steps, t, elapsed);
        }
    }

    double elapsed = (double)(clock() - t_start) / CLOCKS_PER_SEC;
    fprintf(stderr, "Done. %d steps in %.2f s (%.0f steps/s)\n",
            n_steps, elapsed, n_steps / elapsed);

    /* Cleanup */
    free(g_state);
    free(g_k1);
    free(g_k2);
    free(g_k3);
    free(g_k4);
    free(g_tmp);

    return 0;
}
