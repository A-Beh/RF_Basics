% Ali Behfarnia
% Edited 11/2024
% Goal: Type-II PLL, for like Charging pump PLL 

% ===============
% Parameters for Type-II PLL with Active Loop Filter
% ===============
I_p = 10e-6;                         % Charge pump current (A)
R1 = 10e3;                          % Resistance in the filter (ohms)
C1 = 1e-6;                           % Capacitance in the filter (farads)
K_VCO = 2 * pi * 100;         % VCO Gain (rad/s/V) - Tuned for slower dynamics
K_PD = 1;                           % Phase Detector Gain (V/rad)

% ===============
% Simulation Parameters
% ===============
fs = 1e5;                               % Sampling frequency (Hz)
Ts = 1 / fs;                            % Sampling time (s)
T_total = 0.5;                           % Total simulation time (s)
N = round(T_total / Ts);        % Total number of samples
t = (0:N-1) * Ts;                     % Time vector

f_ref_initial = 10;         % Initial frequency of reference signal (Hz)
f_ref_new = 40;             % New frequency of reference signal (Hz)
change_time = T_total/2;         % Time when frequency changes (s)

% ===============
% Reference Phase Signal
% ===============
phi_ref = zeros(1, N);
for n = 2:N
    if t(n) < change_time
        phi_ref(n) = phi_ref(n-1) + 2 * pi * f_ref_initial * Ts;
    else
        phi_ref(n) = phi_ref(n-1) + 2 * pi * f_ref_new * Ts;
    end
end

% ===============
% Filter Coefficients (Discrete-Time)
% ===============
alpha = R1 * C1 / (Ts + R1 * C1);     % RC filter coefficient
beta = I_p / (2 * pi * C1);           % Proportional term

% ===============
% Initialization
% ===============
phi_out = zeros(1, N);      % Output phase of PLL
v_cont = zeros(1, N);       % Control voltage to VCO
v_PD = zeros(1, N);         % Phase detector (voltage)
v_filter = zeros(1, N);     % Filter output (voltage)
phase_error = zeros(1, N);  % Phase error

% ===============
% PLL Simulation
% ===============
for n = 2:N
    % Phase detector: Compute phase error
    phase_error(n) = phi_ref(n) - phi_out(n-1);

    % Phase detector output
    v_PD(n) = K_PD * phase_error(n);

    % Active filter dynamics (discrete implementation)
    v_filter(n) = alpha * v_filter(n-1) + (1 - alpha) * v_PD(n) + beta * phase_error(n);

    % Control voltage to VCO
    v_cont(n) = v_filter(n);

    % Prevent instability by clamping the control voltage
    v_cont(n) = max(min(v_cont(n), 10), -10);

    % VCO: Update output phase
    phi_out(n) = phi_out(n-1) + (K_VCO * v_cont(n) * Ts);
end

% ===============
% Input and Output of PLL
% ===============
ref_signal = sin(phi_ref);  % Reference signal
output_signal = sin(phi_out);  % Output signal

% ===============
% Plot Results
% ===============
figure;
subplot(2, 1, 1);
plot(t, ref_signal, 'k--', 'LineWidth', 4); hold on;
plot(t, output_signal, 'b', 'LineWidth', 1.6);
xlabel('Time (s)', FontSize=16); ylabel('Amplitude', FontSize=16);
title('Reference Signal and PLL Output Signal', FontSize=18);
legend('Reference Signal', 'Output Signal');
grid on;
set(gca, 'FontSize', 16);

subplot(2, 1, 2);
plot(t, phase_error, 'r', 'LineWidth', 2);
xlabel('Time (s)', FontSize=16); ylabel('Phase Error (rad)', FontSize=16);
title('Phase Error over Time', FontSize=18);
grid on;
set(gca, 'FontSize', 16);
