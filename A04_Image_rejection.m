% Edited by A. Behfarnia, 11/2024
%  Hartley Image Rejection Receiver

% ================
% Parameters
% ================
fs = 1e6;                             % Sampling frequency (1 MHz)
T = 1;                                  % Signal duration (1 second)
t = (0:1/fs:T-1/fs);                    % Time vector

% ================
% Signal Parameters
% ================
fc = 200e3;                            % Desired signal center frequency (200 kHz)
f_im = 150e3;                          % Image signal center frequency (150 kHz)
f_LO = 175e3;                          % LO frequency (175 kHz)
BW = 20e3;                              % Signal bandwidth (20 kHz)

% ================
% Generate Desired and Image Signals
% ================
desired_signal = cos(2*pi*fc*t) .* sinc(BW*(t-0.5));  % Desired RF signal with BW
image_signal = cos(2*pi*f_im*t) .* sinc(BW*(t-0.5));  % Image RF signal with BW
RF_signal = desired_signal + image_signal;           % Composite RF signal

% ================
% LO Signals for Quadrature Mixing
% ================
LO_I = cos(2*pi*f_LO*t);      % In-phase LO signal
LO_Q = sin(2*pi*f_LO*t);      % Quadrature LO signal

% ================
% Mix RF signal with LO
% ================
I_signal = RF_signal .* LO_I;  % In-phase component
Q_signal = RF_signal .* LO_Q;  % Quadrature component

% ================
% Apply a 90° Phase Shift to the Q Path
% ================
Q_shifted = imag(hilbert(Q_signal));  % Ideal 90-degree phase shift

% ================
%
% ================

% ================
% Combine I and Q for Image Rejection (Hartley combination)
% ================
combined_signal_1 = I_signal + Q_shifted;  % Desired signal (constructive)
combined_signal_2 = I_signal - Q_shifted;  % Image signal (destructive)

% ================
% Choose the signal that represents the desired signal after rejection
% ================
baseband_signal = combined_signal_1;  % Select the signal after image rejection

% ================
% Frequency Analysis using FFT
% ================
nfft = 2^nextpow2(length(t));
f = linspace(-fs/2, fs/2, nfft);

% ================
% FFT of Signals
% ================
RF_spectrum = fftshift(abs(fft(RF_signal, nfft)));                       % RF signal spectrum
Baseband_spectrum = fftshift(abs(fft(baseband_signal, nfft))); % Baseband spectrum

% ================
% Plots
% ================
figure;

% Before Image Rejection
subplot(2,1,1);
plot(f, 20*log10(RF_spectrum));
title('Spectrum of Composite RF Signal (Before Image Rejection)');
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
grid on;
xlim([-fs/2 fs/2]);

% After Image Rejection
subplot(2,1,2);
plot(f, 20*log10(Baseband_spectrum));
title('Spectrum of Baseband Signal (After Image Rejection)');
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
grid on;
xlim([-50e3 50e3]); % Zoomed around baseband frequencies
