% Ali Behfarnia, Edited 10/2024
% Wi-Fi System Simulation with Phase Noise Impact 
% Modeling of the LO(Local Oscillator) phase noise as a Wiener process (cumulative sum of Gaussian noise)

close all;
clc;
clear all;

% Wi-Fi System Simulation with Phase Noise Impact
% ========================
% Parameters
% ========================
f_c = 2.4e9; 
fs = 4 * f_c;    % Sampling frequency (4 points for the constellation)
ts = 1 / fs;     
t = 0:ts:1.5e-9-ts;   % Time vector, simulate 1.5 nanoseconds

a = cos(2 * pi * f_c * t);   % In-phase (I) component
b = sin(2 * pi * f_c * t);   % Quadrature (Q) component

marker_size = 200; % Marker size for constellation points

% ==================
% Case 1: IQ Imbalance (epsilon = 0.1, theta = pi/6)
% ==================
epsilon1 = 0.2;   
theta1 = 0;    

% Transmit Oscillator without IQ Mismatch
cos_tx = cos(2 * pi * f_c * t);   % Ideal In-phase (I) component
sin_tx = sin(2 * pi * f_c * t);   % Ideal Quadrature (Q) component

% Received signal with IQ imbalance (case 1)
x_bb_i1 = a .* ((1 + epsilon1 / 2) * cos(theta1 / 2)) - b .* ((1 + epsilon1 / 2) * sin(theta1 / 2));
x_bb_q1 = -a .* ((1 - epsilon1 / 2) * sin(theta1 / 2)) + b .* ((1 - epsilon1 / 2) * cos(theta1 / 2));

% ==================
% Case 2: IQ Imbalance (epsilon = 0.2, theta = pi/4)
% ==================
epsilon2 = 0;  
theta2 = pi / 4;   

% Received signal with IQ imbalance (case 2)
x_bb_i2 = a .* ((1 + epsilon2 / 2) * cos(theta2 / 2)) - b .* ((1 + epsilon2 / 2) * sin(theta2 / 2));
x_bb_q2 = -a .* ((1 - epsilon2 / 2) * sin(theta2 / 2)) + b .* ((1 - epsilon2 / 2) * cos(theta2 / 2));

% ==================
% Plotting the Constellations
% ==================
figure;

% Subplot for Case 1
subplot(1, 2, 1);
scatter(cos_tx, sin_tx, marker_size, 'b', 'filled');   % Tx oscillator (ideal)
hold on;
scatter(x_bb_i1, x_bb_q1, marker_size, 'r', 'filled'); % Rx oscillator with IQ mismatch (case 1)
title(['Constellation with \epsilon = ' num2str(epsilon1) ', \theta = ' num2str(theta1)], 'FontSize', 20);
xlabel('In-phase (I)', 'FontSize', 20);
ylabel('Quadrature (Q)', 'FontSize', 20);
ylim([-1.5 1.5]);
grid on;
axis equal;
legend('Tx Oscillator', 'Rx Oscillator with IQ Imbalance', 'FontSize', 14);
hold off;

% Subplot for Case 2
subplot(1, 2, 2);
scatter(cos_tx, sin_tx, marker_size, 'b', 'filled');   % Tx oscillator (ideal)
hold on;
scatter(x_bb_i2, x_bb_q2, marker_size, 'g', 'filled'); % Rx oscillator with IQ mismatch (case 2)
title(['Constellation with \epsilon = ' num2str(epsilon2) ', \theta = ' num2str(theta2)], 'FontSize', 20);
xlabel('In-phase (I)', 'FontSize', 20);
ylabel('Quadrature (Q)', 'FontSize', 20);
ylim([-1.5 1.5]);
grid on;
axis equal;
legend('Tx Oscillator', 'Rx Oscillator with IQ Imbalance', 'FontSize', 14);
hold off;

% ==================
% BER vs SNR Calculation for Four Cases
% ==================

% Simulation parameters for BER
num_symbols = 1e5;          % Number of symbols for BER calculation
snr_values = 0:2:20;        % SNR range in dB

% IQ imbalance parameters for each case
epsilon3 = 0.1; theta3 = pi / 6;   % Case 3: both epsilon and theta non-zero
epsilon_ideal = 0; theta_ideal = 0;% Ideal case: no IQ imbalance

% Initialize BER results
ber_ideal = zeros(1, length(snr_values));
ber_case1 = zeros(1, length(snr_values));
ber_case2 = zeros(1, length(snr_values));
ber_case3 = zeros(1, length(snr_values));

for i = 1:length(snr_values)
    snr = snr_values(i);
    noise_power = 10^(-snr / 10);

    % Generate random BPSK symbols
    a = randi([0, 1], 1, num_symbols) * 2 - 1;
    b = randi([0, 1], 1, num_symbols) * 2 - 1;

    % Ideal case
    x_ideal_i = a;
    x_ideal_q = b;
    
    % Case 1: IQ imbalance with epsilon1 and theta1
    x_bb_i1 = a .* ((1 + epsilon1 / 2) * cos(theta1 / 2)) - b .* ((1 + epsilon1 / 2) * sin(theta1 / 2));
    x_bb_q1 = -a .* ((1 - epsilon1 / 2) * sin(theta1 / 2)) + b .* ((1 - epsilon1 / 2) * cos(theta1 / 2));
    
    % Case 2: IQ imbalance with epsilon2 and theta2
    x_bb_i2 = a .* ((1 + epsilon2 / 2) * cos(theta2 / 2)) - b .* ((1 + epsilon2 / 2) * sin(theta2 / 2));
    x_bb_q2 = -a .* ((1 - epsilon2 / 2) * sin(theta2 / 2)) + b .* ((1 - epsilon2 / 2) * cos(theta2 / 2));
    
    % Case 3: IQ imbalance with epsilon3 and theta3 (both non-zero)
    x_bb_i3 = a .* ((1 + epsilon3 / 2) * cos(theta3 / 2)) - b .* ((1 + epsilon3 / 2) * sin(theta3 / 2));
    x_bb_q3 = -a .* ((1 - epsilon3 / 2) * sin(theta3 / 2)) + b .* ((1 - epsilon3 / 2) * cos(theta3 / 2));
    
    % Add noise
    noise_i = sqrt(noise_power / 2) * randn(1, num_symbols);
    noise_q = sqrt(noise_power / 2) * randn(1, num_symbols);

    % Received symbols for each case
    r_ideal_i = x_ideal_i + noise_i;
    r_ideal_q = x_ideal_q + noise_q;
    r_i1 = x_bb_i1 + noise_i;
    r_q1 = x_bb_q1 + noise_q;
    r_i2 = x_bb_i2 + noise_i;
    r_q2 = x_bb_q2 + noise_q;
    r_i3 = x_bb_i3 + noise_i;
    r_q3 = x_bb_q3 + noise_q;

    % Detection and BER calculation
    ber_ideal(i) = mean(sign(r_ideal_i) ~= a) + mean(sign(r_ideal_q) ~= b) / 2;
    ber_case1(i) = mean(sign(r_i1) ~= a) + mean(sign(r_q1) ~= b) / 2;
    ber_case2(i) = mean(sign(r_i2) ~= a) + mean(sign(r_q2) ~= b) / 2;
    ber_case3(i) = mean(sign(r_i3) ~= a) + mean(sign(r_q3) ~= b) / 2;
end

% Plotting BER vs SNR
figure;
semilogy(snr_values, ber_ideal, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Ideal (No IQ Imbalance)');
hold on;
semilogy(snr_values, ber_case1, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', ['\epsilon = ' num2str(epsilon1) ', \theta = ' num2str(theta1)]);
semilogy(snr_values, ber_case2, '-^', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', ['\epsilon = ' num2str(epsilon2) ', \theta = ' num2str(theta2)]);
semilogy(snr_values, ber_case3, '-d', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', ['\epsilon = ' num2str(epsilon3) ', \theta = ' num2str(theta3)]);
grid on;
xlabel('SNR (dB)', 'FontSize', 14);
ylabel('Bit Error Rate (BER)', 'FontSize', 14);
title('BER vs SNR for Various IQ Imbalance Conditions', 'FontSize', 16);
legend('Location', 'southwest', 'FontSize', 16);
hold off;
