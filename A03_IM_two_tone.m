% Ali Behfarnia
% Edited 2024
% Goal: two-tone test IP3 calculation

% ===========
% Parameters
% ===========
Fs = 1e6;                  % Sampling frequency (1 MHz)
T = 1;                        % Signal duration (1 second)
t = 0:1/Fs:T-1/Fs;      % Time vector

f_1 = 2.4e5;               % Frequency of first tone (240 kHz, scaled for demonstration)
f_2 = 2.401e5;           % Frequency of second tone, slightly offset (1 kHz apart)

Pwr_tone = -10;                                      % Tone power (dBm)
Amp_tone = 10^((Pwr_tone - 30) / 20);  % Convert dBm to linear scale (voltage)

% ===========
% Generate two-tone signal
% ===========
tone1 = Amp_tone * cos(2*pi*f_1*t);
tone2 = Amp_tone * cos(2*pi*f_2*t);
input_signal = tone1 + tone2;

% ===========
%  Coefficients for main and the 3rd-order non-linearity (cubic non-linear model)
% ===========
a1 = 1;            % Linear gain
a3 = 1e-5;       % 3rd-order non-linearity coefficient
output_signal = a1 * input_signal + a3 * input_signal.^3;

% ===========
% Perform FFT to observe frequency spectrum
% ===========
N = length(t);
f = (-N/2:N/2-1)*(Fs/N); 
Input_FFT = abs(fftshift(fft(input_signal)))/N;
Output_FFT = abs(fftshift(fft(output_signal)))/N;

% ===========
% Plot: frequency spectrum
% ===========
figure;
subplot(2,1,1);
plot(f, 20*log10(Input_FFT));
title('Frequency Spectrum of Input Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

subplot(2,1,2);
plot(f, 20*log10(Output_FFT));
title('Frequency Spectrum of Output Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

% ===========
% Calculating Third-Order Intercept Point (IP3)
% ===========
% Convert frequencies to corresponding indices in FFT array
f1_index = round(N/2 + f_1*(N/Fs));
f2_index = round(N/2 + f_2*(N/Fs));

% IMD (Inter-Modulation Distortion) frequencies
IMD_freq1_index = round(N/2 + (2*f_1 - f_2)*(N/Fs));
IMD_freq2_index = round(N/2 + (2*f_2 - f_1)*(N/Fs));

% Power at fundamentals and IMD products
P_f1 = Output_FFT(f1_index); % Fundamental at f1
P_f2 = Output_FFT(f2_index); % Fundamental at f2
P_IMD1 = Output_FFT(IMD_freq1_index); % IMD product at 2f1-f2
P_IMD2 = Output_FFT(IMD_freq2_index); % IMD product at 2f2-f1

P_input = Pwr_tone;  
IP3_dB = P_input + (20*log10(P_f1) - 20*log10(P_IMD1)); 
fprintf('Third-Order Intercept Point (IP3): %.2f dBm \n', IP3_dB);

