% Ali Behfarnia
% Edited 2024
% Goal: Dual Downconversion Heterodyne System with Mixing Spurs
%           This includes 2 interferers and 2nd harmonic of LO2

close all;
clear all;
clc;

% ===========
% Parameters
% ===========
Fs = 10e9;                            % Sampling frequency (10 GHz for high-resolution simulation)
t = 0:1/Fs:1e-6;                     % Time: 1 micro

% ===========
% Signal and interference parameters
% ===========
f_signal = 2.4e9;                   % Desired signal center frequency (2.4 GHz)
signal_BW = 20e6;             % Bandwidth of the desired signal (200 MHz)
f_interferer1 = 2.76e9;          % f_c for Interferer 1
f_interferer2 = 2.8e9;            % f_c for Interferer 2
interferer_BW = 10e6;          % BW of the interferers (10 MHz)

% ===========
% LO frequencies
% ===========
LO1_freq = 1.98e9; % First Local Oscillator frequency
LO2_freq = 400e6; % Second Local Oscillator frequency (LO2)
LO2_harmonic_freq = 2 * LO2_freq; % Second harmonic of LO2 (800 MHz)

% ===========
% Generate signals
% ===========
desired_signal = cos(2*pi*f_signal*t) .* sinc(signal_BW*(t - mean(t))); 
interferer_1 =  cos(2*pi*f_interferer1*t) .* sinc(interferer_BW*(t - mean(t))); 
interferer_2 =  cos(2*pi*f_interferer2*t) .* sinc(interferer_BW*(t - mean(t))); 

received_signal = desired_signal + interferer_1 + interferer_2;

% ===========
% Plot Received Spectrum (With Interferers)
% ===========
figure;
n_FFT = 2^nextpow2(length(received_signal)); % FFT points
f = (-n_FFT/2:n_FFT/2-1)*(Fs/n_FFT);               % Frequency axis
received_spectrum = fftshift(abs(fft(received_signal, n_FFT)));
plot(f/1e9, received_spectrum, 'LineWidth',2); 
title('Received Spectrum (With Interferers)');
xlabel('Frequency (GHz)');
ylabel('Amplitude');
xlim([2 3]);
grid on;


% ===========
% Step 1: Amplify using LNA
% ===========
LNA_gain = 20;                      
amplified_signal = received_signal * 10^(LNA_gain/20);

% ===========
% Step 2: First Downconversion (Mix with LO1)
% ===========
LO1 = cos(2*pi*LO1_freq*t);     % LO1 signal
IF1 = amplified_signal .* LO1; % First mixing with interferers

% ===========
% Step 3: Bandpass Filter at IF1 (420 MHz)
% ===========
BPF_bandwidth = 800e6; % BW of BPF (it does pass the interferers) 
[b, a] = butter(4, [420e6-BPF_bandwidth/2, 420e6+BPF_bandwidth/2]/(Fs/2));
filtered_IF1 = filter(b, a, IF1);

% ===========
% Step 4: Second Downconversion (Mix with LO2 and its 2nd harmonic)
% ===========
LO2 = cos(2*pi*LO2_freq*t); % LO2 signal
LO2_harmonic = cos(2*pi*LO2_harmonic_freq*t); % 2nd harmonic of LO2 signal

% Case 1: With 2nd Harmonic
baseband_signal_with_harmonic = filtered_IF1 .* LO2 + filtered_IF1 .* LO2_harmonic; % Second mixing with interferers and its harmonic

% Case 2: Ideal Case (No 2nd Harmonic)
baseband_signal_ideal = filtered_IF1 .* LO2; % Second mixing only with LO2 (no harmonic)

% ===========
% Step 5: Low-Pass Filter (LPF) after second mixer to remove high frequencies
% ===========
LPF_cutoff = 40e6;                                               % LPF cutoff frequency
[b_lpf, a_lpf] = butter(4, LPF_cutoff/(Fs/2));          
filtered_output_with_harmonic = filter(b_lpf, a_lpf, baseband_signal_with_harmonic); 
filtered_output_ideal = filter(b_lpf, a_lpf, baseband_signal_ideal); 

% ===========
% Plots into Subplots (2x2)
% ===========
figure;

% Define font size parameters
axis_font_size = 14;
label_font_size = 16;

% Subplot 1: Time Domain Signal (With Harmonic)
subplot(2, 2, 2);
plot(t*1e6, filtered_output_with_harmonic, 'r', 'LineWidth', 2); % Time in microseconds
title('Time Domain (With LO2 2nd Harmonic)', 'FontSize', label_font_size);
xlabel('Time (µs)', 'FontSize', label_font_size);
ylabel('Amplitude', 'FontSize', label_font_size);
grid on;
set(gca, 'FontSize', axis_font_size);

% Subplot 2: Time Domain Signal (Ideal Case)
subplot(2, 2, 1);
plot(t*1e6, filtered_output_ideal, 'LineWidth', 2); 
title('Time Domain (Ideal, No LO2 2nd Harmonic)', 'FontSize', label_font_size);
xlabel('Time (µs)', 'FontSize', label_font_size);
ylabel('Amplitude', 'FontSize', label_font_size);
grid on;
set(gca, 'FontSize', axis_font_size);

% Subplot 3: Frequency Domain Spectrum (With Harmonic)
subplot(2, 2, 4);
output_spectrum_with_harmonic = fftshift(abs(fft(filtered_output_with_harmonic, n_FFT)));
plot(f/1e9, output_spectrum_with_harmonic, 'r', 'LineWidth', 2); % Convert frequency to GHz
title('Frequency Domain (With LO2 2nd Harmonic)', 'FontSize', label_font_size);
xlabel('Frequency (GHz)', 'FontSize', label_font_size);
ylabel('Amplitude', 'FontSize', label_font_size);
xlim([-0.1 0.1]); 
grid on;
set(gca, 'FontSize', axis_font_size);

% Subplot 4: Frequency Domain Spectrum (Ideal Case)
subplot(2, 2, 3);
output_spectrum_ideal = fftshift(abs(fft(filtered_output_ideal, n_FFT)));
plot(f/1e9, output_spectrum_ideal, 'LineWidth', 2); 
title('Frequency Domain (Ideal, No LO2 2nd Harmonic)', 'FontSize', label_font_size);
xlabel('Frequency (GHz)', 'FontSize', label_font_size);
ylabel('Amplitude', 'FontSize', label_font_size);
xlim([-0.1 0.1]); 
grid on;
set(gca, 'FontSize', axis_font_size);

