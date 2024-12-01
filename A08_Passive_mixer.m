% Ali Behfarnia
% Edited 11/2024
% Goal: Single/Double-balanced passive mixers
% Effect of the LO Feedthrough 

clc;
close all;

% ===========
% Parameters
% ===========
fs = 1e6;                                                           % Sampling frequency (Hz)
t = 0:1/fs:1e-3;                                                  % Time vector
f_RF = 100e3;                                                  % RF signal frequency (Hz)
f_LO = 150e3;                                                  % LO signal frequency (Hz)
A_RF = 1;                                                         % RF signal amplitude
A_LO = 1;                                                         % LO signal amplitude
LO_leakage_ratios = [0, 0.2, 0.5, 1, 2];            % Extended LO leakage ratios

% ===========
% Signals
% ===========
RF = A_RF * cos(2 * pi * f_RF * t);                    % RF signal
LO = A_LO * square(2 * pi * f_LO * t);              % LO signal (square wave)

% ===========
% Double-Balanced Mixer
% ===========
LO_pos = max(LO, 0);                        % Positive LO polarity
LO_neg = max(-LO, 0);                       % Negative LO polarity
double_balanced_output = (RF .* LO_pos) - (RF .* LO_neg);

% ===========
% FFT for double-balanced mixer
% ===========
N = length(t);
frequencies = linspace(-fs/2, fs/2, N);
double_fft = fftshift(abs(fft(double_balanced_output)) / N);

% ===========
% Plot: Double_balanced mixer
% ===========
figure;
plot(frequencies/1e3, double_fft, 'k', 'LineWidth', 2, 'DisplayName', 'Double-Balanced (No Leakage)'); 
hold on;

% ===========
% Single-balanced mixer with varying LO leakage ratios
% ===========
colors = lines(length(LO_leakage_ratios));                 
single_ffts = zeros(length(LO_leakage_ratios), N);   
for idx = 1:length(LO_leakage_ratios)
    ratio = LO_leakage_ratios(idx);
   
    single_balanced_raw = RF .* LO;                  % Single-balanced Mixer Output
    LO_leakage = ratio * LO;                                % Simulate LO leakage
    single_balanced_output = single_balanced_raw + LO_leakage; % Total
    
    % FFT for single-balanced mixer
    single_fft = fftshift(abs(fft(single_balanced_output)) / N);
    single_ffts(idx, :) = single_fft;        % Store for later zoomed-in plot
    
    % Plot the spectrum
    plot(frequencies/1e3, single_fft, 'Color', colors(idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Single-Balanced (Leakage %.1f)', ratio));
end

% Plot settings
title('Frequency-Domain Comparison: Double vs. Single Balanced Mixer');
xlabel('Frequency (kHz)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
xlim([-300 300])
legend('show', 'Location', 'best');
grid on;
set(gca, 'FontSize', 20);


% ===========
% Plot: Zoomed-in
% ===========
figure;
hold on;
LO_peak_range = [f_LO - 50, f_LO + 50]; % Extended frequency range around the LO peak
for idx = 1:length(LO_leakage_ratios)
    plot(frequencies/1e3, single_ffts(idx, :), 'Color', colors(idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Single-Balanced (Leakage %.1f)', LO_leakage_ratios(idx)));
end

plot(frequencies/1e3, double_fft, 'k', 'LineWidth', 2, 'DisplayName', 'Double-Balanced (No Leakage)');
xlim(LO_peak_range / 1e3);
title('Zoomed-In: LO Frequency Region with Extended Bandwidth');
xlabel('Frequency (kHz)' , 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
legend('show', 'Location', 'best');
grid on;
set(gca, 'FontSize', 20);
