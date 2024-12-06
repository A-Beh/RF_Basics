% Ali Behfarnia
% Edited 11/2024
% Goal: Type I PLL

% Notes:
% 1) Two Ref. Frequencies 50 and 100 are followed by the output
% 2) An RC Low pass Filter is used fo loop filter

% ===============
% Parameters for Type-I PLL with RC filter
% ===============
K_PD = 1;                % Phase Detector Gain (V/rad)
R1 = 1;                     % Resistance in the RC filter (ohms)
C1 = 1e-3;                % Capacitance in the RC filter (farads)
K_VCO = 200;          % VCO Gain (rad/s/V)

% ===============
% Simulation parameters
%===============
fs = 1e4;                     % Sampling frequency (Hz)
Ts = 1/fs;                    % Sampling time (s)
N = 1000;                   % Number of samples
t = (0:N-1) * Ts;           % Time vector

% ===============
% Signals
% ===============
f_ref_initial = 50;         % Initial frequency of reference signal (Hz)
f_ref_new = 100;         % New frequency of reference signal (Hz)
change_time = 0.05;   % Time when the frequency changes

% ===============
%  Phase of Signal
% ===============
phi_ref = zeros(1, N);
for n = 2:N
    if t(n) < change_time
        phi_ref(n) = phi_ref(n-1) + 2*pi*f_ref_initial*Ts;      % Initial frequency
    else
        phi_ref(n) = phi_ref(n-1) + 2*pi*f_ref_new*Ts;        % New frequency
    end
end

% ===============
%  Initialization
% ===============
phi_out = zeros(1, N);               % Output phase of PLL
v_cont = zeros(1, N);                % Control voltage fed to the VCO.
v_PD = zeros(1, N);                  % Phase detector (voltage)
v_filter = zeros(1, N);                % Voltage across the RC = v_cont
phase_error = zeros(1, N);       % Phase error

% ===============
%  RC filter coefficient (alpha)
% ===============
alpha = R1 * C1 / (Ts + R1 * C1);

% ===============
% PLL simulation
% ===============
for n = 2:N
    % Phase detector: Compute phase error
    phase_error(n) = phi_ref(n) - phi_out(n-1);

     % Phase detector output
    v_PD(n) = K_PD * phase_error(n);

    % RC filter dynamics
    v_filter(n) = alpha * v_filter(n-1) + (1 - alpha) * v_PD(n);
    
    % Control voltage to VCO
    v_cont(n) = v_filter(n); 

    % VCO: Update output phase
    phi_out(n) = phi_out(n-1) + (K_VCO * v_cont(n) * Ts);
end

% ===============
% Input and Output of PLL
% ===============
ref_signal = sin(phi_ref); % Reference signal
output_signal = sin(phi_out); % Output signal

% ===============
% Plot results
% ===============
figure;
subplot(2, 1, 1);
plot(t, ref_signal, 'k', 'LineWidth', 1.2); hold on;
plot(t, output_signal, 'b', 'LineWidth', 1.2);
xlabel('Time (s)', FontSize=16); ylabel('Amplitude' , FontSize=16);
title('Reference Signal and PLL Output Signal', FontSize=18);
legend('Reference Signal', 'Output Signal');
grid on;
set(gca, 'FontSize', 16);

subplot(2, 1, 2);
plot(t, phase_error, 'r', 'LineWidth', 1.4);
xlabel('Time (s)', FontSize=16); ylabel('Phase Error (rad)',FontSize=16);
title('Phase Error over Time', FontSize=18);
grid on;
set(gca, 'FontSize', 16);