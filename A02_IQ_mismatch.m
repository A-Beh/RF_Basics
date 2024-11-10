% Editted By Ali Behfarnia
% Editted 10/2024
% Goal: I/Q Imbalance Representation

clc; clear all; close all;

% Sampling parameters for Wi-Fi-like system
f_c = 2.4e9; 
fs = 8 * f_c;                                                                               % Sampling frequency, e.g., 10 GHz
ts = 1/fs;                                                                                   % Sampling time
t = 0:ts:2e-9-ts;                                                                         % Time vector, simulate 2 nanoseconds


% ==================
% Step 1: Perfect  Transmit Oscillator
% ==================

% Transmit oscillator (perfect Tx oscillator without I/Q mismatch)
cos_tx = cos(2*pi*f_c*t);                                                          % In-phase (I) component
sin_tx = sin(2*pi*f_c*t);                                                            % Quadrature (Q) component

% Plotting the Lissajous figure for the Tx oscillator (constellation)
% figure;
% scatter(cos_tx, sin_tx, 'b', 'filled');                                            % Tx oscillator (circle)
% title('Transmit Oscillator (Constellation)');
% xlabel('In-phase (I)');
% ylabel('Quadrature (Q)');
% grid on;
% axis equal;                                                                     

% ==================
% Step 2: Received Signal with IQ Imbalance
% ==================

% Parameters for IQ imbalance
psi = pi/5;                                                                             % Phase offset (radians) for Q
alpha = 0.9;                                                                         % Amplitude imbalance factor for I 

% Define the received carrier oscillator model
cos_rx = alpha * cos(2 * pi * f_c * t);                                    % Scaled I component
sin_rx = sin(2 * pi * f_c * t + psi);                                         % Q component with phase offset


% Plotting the Lissajous figure for the Rx oscillator (constellation)
markerSize = 100;
figure;
scatter(cos_rx, sin_rx, markerSize, 'r', 'filled');                                            % Rx oscillator (ellipsoid)
title('Received Oscillator (Constellation)',  'FontSize', 22);
xlabel('In-phase (I)', 'FontSize', 18);
ylabel('Quadrature (Q)', 'FontSize', 18);
grid on;
axis equal;                                                                               % Ensure equal scaling on both axes
hold on;
scatter(cos_tx, sin_tx, markerSize, 'b', 'filled'); % Overlay the Tx constellation
legend('Rx Oscillator', 'Tx Oscillator', 'FontSize', 18);
hold off;