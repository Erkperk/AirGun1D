#!/usr/bin/env python3
"""Plot airgun bubble oscillation results and chamber pressure space-time."""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_bubble(fname):
    data = np.loadtxt(fname, comments='#')
    t = data[:, 0]
    R = data[:, 1]
    Rdot = data[:, 2]
    p_bub = data[:, 3]
    p_ac = data[:, 4]

    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t * 1e3, R, 'b-', linewidth=1.5)
    axes[0].set_ylabel('R [m]')
    axes[0].set_title('Bubble oscillations')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t * 1e3, Rdot, 'b-', linewidth=1.5)
    axes[1].set_ylabel(r'$\dot{R}$ [m/s]')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t * 1e3, p_bub / 1e5, 'b-', linewidth=1.5)
    axes[2].set_ylabel('Bubble pressure [bar]')
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(t * 1e3, p_ac / 1e5, 'b-', linewidth=1.5)
    axes[3].set_ylabel('Acoustic pressure [bar]')
    axes[3].set_xlabel('Time [ms]')
    axes[3].grid(True, alpha=0.3)

    plt.tight_layout()

    out = fname.rsplit('.', 1)[0] + '.png'
    plt.savefig(out, dpi=150)
    print(f'Saved: {out}')
    plt.close()

def plot_chamber(chamber_fname):
    from matplotlib.colors import LogNorm

    with open(chamber_fname, 'r') as f:
        header = f.readline()

    # Parse grid coordinates from header: "# x: x0 x1 x2 ..."
    x = np.array([float(v) for v in header.split(':')[1].split()])

    # Load space-time data: each row is [t, p0, p1, ...]
    data = np.loadtxt(chamber_fname, comments='#')
    t = data[:, 0]
    P = data[:, 1:]  # pressure field (n_times x n_x)
    P_bar = P / 1e5

    # Trim trailing frozen state (after cutoff, RHS=0 so pressure is constant)
    dp = np.max(np.abs(np.diff(P_bar, axis=0)), axis=1)
    active = np.where(dp > 1e-6)[0]
    if len(active) > 0:
        i_end = min(active[-1] + 20, len(t))
    else:
        i_end = len(t)

    t_plot = t[:i_end]
    P_plot = P_bar[:i_end]

    # Clamp for log scale
    P_plot = np.clip(P_plot, 0.1, None)

    fig, ax = plt.subplots(figsize=(10, 6))
    pcm = ax.pcolormesh(x, t_plot * 1e3, P_plot, shading='auto',
                        cmap='inferno', norm=LogNorm(vmin=0.1, vmax=P_plot.max()))
    cb = fig.colorbar(pcm, ax=ax)
    cb.set_label('Pressure [bar]')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('Time [ms]')
    ax.set_title('Chamber pressure (space-time)')

    plt.tight_layout()

    out = chamber_fname.rsplit('.', 1)[0] + '.png'
    plt.savefig(out, dpi=150)
    print(f'Saved: {out}')
    plt.close()

def main():
    if len(sys.argv) > 1:
        fname = sys.argv[1]
    else:
        fname = 'airgun_output.dat'

    plot_bubble(fname)

    # Look for chamber file
    base = fname.rsplit('.', 1)[0]
    chamber_fname = base + '_chamber.dat'
    if os.path.exists(chamber_fname):
        plot_chamber(chamber_fname)

if __name__ == '__main__':
    main()
