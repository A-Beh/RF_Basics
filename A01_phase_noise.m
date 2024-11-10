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
M = 4;                                             % M-QAM Modulation, M = 4
N = 1000;                                       % Number of symbols
Fs = 40e6;                                      % Sampling frequency (Wi-Fi channel bandwidth = 20 MHz)
f_c = 2.4e9;                                    % Carrier frequency (2.4 GHz for Wi-Fi)
SNR_dB = 25;                                % Signal-to-Noise Ratio in dB
Ts = 1/Fs;                                       % Symbol period

% Phase noise variance Vector
phase_noise_variances = [1e-6 1e-5 01e-4 1e-3];

% ========================
% QAM Modulation
% ========================
data = randi([0 M-1], N, 1);                 % Random data symbols
tx_symbols = qammod(data, M, 'UnitAveragePower', true);  % M-QAM modulation

% ========================
% Plot Constellation Diagrams (With LO Phase Noise)
% ========================
figure;

for i = 1:length(phase_noise_variances)
    phase_noise_variance = phase_noise_variances(i);
    
    % LO phase noise
    delta_phase = sqrt(phase_noise_variance) * randn(N, 1);  % Gaussian increments
    phase_noise = cumsum(delta_phase);  % Cumulative sum gives a Wiener process (random walk)

    % Phase noise added
    tx_symbols_with_phase_noise = tx_symbols .* exp(1j * phase_noise);  % Apply cumulative phase noise

    % AWGN
    rx_symbols_noiseless = tx_symbols_with_phase_noise; 
    rx_symbols_noisy = awgn(rx_symbols_noiseless, SNR_dB, 'measured'); 

    scatterplot = @(x) scatter(real(x), imag(x));

    % Plot using a 2*2 
    subplot(2, 2, i);  
    scatterplot(rx_symbols_noisy);
    xlabel('In-phase')
    ylabel('Quadrature')
    title(['4-QAM Constellation with LO Phase Noise (Variance = ', num2str(phase_noise_variance), ')']);
    grid on;
end
