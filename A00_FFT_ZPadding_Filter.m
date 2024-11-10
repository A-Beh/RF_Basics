% Ali Behfarnia
% Created 2017, Edited 2024
% Goal: Sample FFT Plot with Zero Padding and Filtering
close all;
clear all;
clc;

% ===========
% Parameters
% ===========
f_sample = 500;                       % Sampling Frequency (Hz)
t_sample = 1/f_sample;                % Sampling time (s)
duration = 0.5;                             % Signal duration (s)
N = f_sample * duration;              % Total number of samples
time = 0:t_sample:duration - t_sample; % Time Vector

% ===========
% Defining two signals
% ===========
a_1 = 1; f_1 = 100; phi_1 = 0.6;
S_1 = a_1*cos(2*pi*f_1*time + phi_1);
a_2 = 2; f_2 = 150; phi_2 = -0.8;
S_2 = a_2*cos(2*pi*f_2*time + phi_2);
S = S_1 + S_2;  % Combined signal

% ===========
% Zero Padding
% ===========
zero_padding_factor = 2;             % You can adjust this factor to add more padding
N_padded = N * zero_padding_factor;  % New length with zero padding
S_padded = [S, zeros(1, N_padded - N)];  % Signal with zero padding

% ===========
% FFT - Without Zero Padding
% ===========
S_freq = fft(S);
S_oneside = S_freq(1:N/2);
f = f_sample * (0:(N/2-1)) / N;
S_mag = abs(S_oneside) / (N/2);

% ===========
% FFT - With Zero Padding
% ===========
S_freq_padded = fft(S_padded);
S_oneside_padded = S_freq_padded(1:N_padded/2);
f_padded = f_sample * (0:(N_padded/2-1)) / N_padded;
S_mag_padded = abs(S_oneside_padded) / (N/2);

% ===========
% Digital Filters
% ===========
% Design a "very good" filter (high order, sharp cut-off)
order_good = 10;
cutoff_freq_good = 120; % Cutoff just above 100 Hz
[b_good, a_good] = butter(order_good, cutoff_freq_good/(f_sample/2), 'low');

% Design a "not-so-good" filter (low order, less sharp)
order_bad = 3;
cutoff_freq_bad = 120; % Same cutoff but lower order
[b_bad, a_bad] = butter(order_bad, cutoff_freq_bad/(f_sample/2), 'low');

% Apply the filters
S_filtered_good = filter(b_good, a_good, S); % Filtered by good filter
S_filtered_bad = filter(b_bad, a_bad, S);    % Filtered by bad filter

% FFT of filtered signals
S_freq_filtered_good = fft(S_filtered_good);
S_oneside_filtered_good = S_freq_filtered_good(1:N/2);
S_mag_filtered_good = abs(S_oneside_filtered_good) / (N/2);

S_freq_filtered_bad = fft(S_filtered_bad);
S_oneside_filtered_bad = S_freq_filtered_bad(1:N/2);
S_mag_filtered_bad = abs(S_oneside_filtered_bad) / (N/2);

% ===========
% Plotting Time Domain and Frequency Domain
% ===========
figure;
sgtitle('Comparison of Zero Padding, Good vs. Bad Filters in Time and Frequency Domains', 'FontSize', 16);

% Plot Time-domain Signal (Original)
subplot(5, 2, 1);
plot(time, S);
xlabel('Time (s)', 'FontSize', 12);
xlim([0 duration]);  
ylabel('Amplitude', 'FontSize', 12);
title('Time-domain Plot (Original)', 'FontSize', 14);
grid on;

% Plot Frequency-domain (Without Zero Padding)
subplot(5, 2, 2);
plot(f, S_mag);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Frequency-domain Plot (No Zero Padding)', 'FontSize', 14);
grid on;

% Plot Time-domain Signal with Zero Padding
subplot(5, 2, 3);
time_padded = 0:t_sample:(length(S_padded)-1)*t_sample; % New time vector for padded signal
plot(time_padded, S_padded);
xlabel('Time (s)', 'FontSize', 12);
xlim([0 duration*zero_padding_factor]);  
ylabel('Amplitude', 'FontSize', 12);
title('Time-domain Plot (With Zero Padding)', 'FontSize', 14);
grid on;

% Plot Frequency-domain (With Zero Padding)
subplot(5, 2, 4);
plot(f_padded, S_mag_padded);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Frequency-domain Plot (With Zero Padding)', 'FontSize', 14);
grid on;

% Plot Bode (Magnitude Response) for "Good" Filter
subplot(5, 2, 5);
[H_good, f_response_good] = freqz(b_good, a_good, N, f_sample);
plot(f_response_good, abs(H_good));
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Magnitude Response (Good Filter)', 'FontSize', 14);
grid on;

% Plot Bode (Magnitude Response) for "Not-So-Good" Filter
subplot(5, 2, 6);
[H_bad, f_response_bad] = freqz(b_bad, a_bad, N, f_sample);
plot(f_response_bad, abs(H_bad));
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Magnitude Response (Bad Filter)', 'FontSize', 14);
grid on;

% Plot Time-domain Signal After "Good" Filter
subplot(5, 2, 7);
plot(time, S_filtered_good);
xlabel('Time (s)', 'FontSize', 12);
xlim([0 duration]);  
ylabel('Amplitude', 'FontSize', 12);
title('Time-domain Plot (After Good Filter)', 'FontSize', 14);
grid on;

% Plot Frequency-domain Signal After "Good" Filter
subplot(5, 2, 8);
plot(f, S_mag_filtered_good);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Frequency-domain Plot (After Good Filter)', 'FontSize', 14);
grid on;

% Plot Time-domain Signal After "Not-So-Good" Filter
subplot(5, 2, 9);
plot(time, S_filtered_bad);
xlabel('Time (s)', 'FontSize', 12);
xlim([0 duration]);  
ylabel('Amplitude', 'FontSize', 12);
title('Time-domain Plot (After Bad Filter)', 'FontSize', 14);
grid on;

% Plot Frequency-domain Signal After "Not-So-Good" Filter
subplot(5, 2, 10);
plot(f, S_mag_filtered_bad);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Frequency-domain Plot (After Bad Filter)', 'FontSize', 14);
grid on;
