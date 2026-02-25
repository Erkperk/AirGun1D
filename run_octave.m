%% AirGun1D - Octave headless run script
% Equivalent to example.m but without GUI plotting

addpath('SBPSAT/');
addpath('SeismicAirgunCode/');
addpath('sbplib/');

% define discretization
nx = 10; % number of grid points per 1m of source length

% source properties (SI units)
src_pressure = 164e5;   % source pressure [Pa] (164 bar)
src_length = 10;        % source length [m]
src_area = 1.0;         % cross-sectional area [m^2]
src_depth = 7.5;        % source depth [m]

% run simulation
fprintf('Running AirGun1D simulation (nx=%d, %.1f bar, %.1f m, %.2f m^2, %.1f m depth)...\n', ...
    nx, src_pressure*1e-5, src_length, src_area, src_depth);
tic
sol = runEulerCode(nx, src_pressure, src_length, src_area, src_depth);
elapsed = toc;
fprintf('Simulation completed in %.2f seconds.\n', elapsed);

% extract results
t = sol.x;
R = sol.y(1,:);   % bubble radius (m)
U = sol.y(2,:);   % bubble wall velocity (m/s)

fprintf('\n--- Bubble Dynamics ---\n');
fprintf('Time steps:       %d\n', length(t));
fprintf('Time range:       [%.4f, %.4f] s\n', t(1), t(end));
fprintf('Max bubble radius: %.4f m\n', max(R));
fprintf('Min bubble radius: %.4f m\n', min(R(R > 0)));

% compute acoustic pressure
r = 75;          % distance from source to receiver (m)
c_inf = 1482;    % speed of sound in water (m/s)
rho_inf = 1000;  % density of water (kg/m^3)

[~, solDY] = deval(sol, t);
A = solDY(2,:);  % bubble wall acceleration (m/s^2)

[tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r);
pDirBarM = pDir * 1e-5 * r;  % convert to bar*m

fprintf('\n--- Acoustic Pressure (direct arrival at %d m) ---\n', r);
fprintf('Peak pressure:    %.4f bar*m\n', max(pDirBarM));
fprintf('Min pressure:     %.4f bar*m\n', min(pDirBarM));

% extract chamber pressure
x = [0:ceil(src_length*nx)]./nx;
rho_ch = zeros(length(x), length(t));
rhov_ch = zeros(length(x), length(t));
e_ch = zeros(length(x), length(t));

for i = 1:length(x)
    rho_ch(i,:) = sol.y(3*i+2,:);
    rhov_ch(i,:) = sol.y(3*i+3,:);
    e_ch(i,:) = sol.y(3*i+4,:);
end

gamma = 1.4;
v_ch = rhov_ch ./ rho_ch;
p_ch = (gamma-1) * (e_ch - 0.5*rho_ch.*v_ch.^2);

pa2bar = 1e-5;
fprintf('\n--- Chamber Pressure ---\n');
fprintf('Initial pressure: %.1f bar\n', p_ch(1,1)*pa2bar);
fprintf('Final pressure:   %.1f bar\n', mean(p_ch(:,end))*pa2bar);

% determine when airgun actually stops firing
v_port = rhov_ch(end,:) ./ rho_ch(end,:);
idx_stop = find(v_port < 0, 1);
if ~isempty(idx_stop)
    fprintf('\n--- Firing Duration ---\n');
    fprintf('Airgun stops at:  %.4f s (flow reversal to inflow)\n', t(idx_stop));
else
    fprintf('\n--- Firing Duration ---\n');
    fprintf('Airgun fires for full cutoff window (%.2f s)\n', 0.30);
end

fprintf('\nSimulation successful!\n');

%% Generate PNG plots

% Plot 1: Bubble radius vs time
fig1 = figure('visible', 'off');
plot(t, R, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Bubble radius (m)');
title('Bubble Radius vs Time');
grid on;
print(fig1, 'bubble_radius.png', '-dpng', '-r150');
fprintf('Saved bubble_radius.png\n');

% Plot 2: Acoustic pressure vs time
fig2 = figure('visible', 'off');
plot(tDir, pDirBarM, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Pressure (bar\cdotm)');
title('Acoustic Pressure vs Time');
grid on;
print(fig2, 'acoustic_pressure.png', '-dpng', '-r150');
fprintf('Saved acoustic_pressure.png\n');

% Plot 3: Chamber pressure space-time
fig3 = figure('visible', 'off');
[T, X] = meshgrid(t, x);
h = surf(X, T*1000, p_ch*pa2bar);
view(2); shading interp;
xlabel('Position (m)');
ylabel('Time (ms)');
cb = colorbar; ylabel(cb, 'bar');
title('Chamber Pressure (space-time)');
ylim([0 20]);
print(fig3, 'chamber_pressure.png', '-dpng', '-r150');
fprintf('Saved chamber_pressure.png\n');
