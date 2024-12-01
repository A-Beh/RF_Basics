% Ali Behfarnia
% Edited 11/2024
% Goal: A cascode of CG LNA, output frequency response

% ===========
% Parameters
% ===========
gm2 = 10e-3;                       % Transconductance of M2 (S)
ro1 = 20e3;                          % Output resistance of M1 (ohms)
Cx = 1e-12;                          % Parasitic capacitance (F)
R1 = 1e3;                             % Load resistance (ohms)
Rs = 50;                                % Antenna resistor

frequencies = logspace(6, 12, 1000);   % From 1 MHz to 1 THz
omega = 2 * pi * frequencies;               % Angular frequency (rad/s)
s = 1j * omega;                                     % Complex frequency variable

% ===========
% Transfer function, abs(Vout/Vn2)
% ===========
Vout = 2 * ro1 * Cx * s + 1;                                                    % Numerator
Vn2 = 2 * ro1 * Cx * s + 2 * gm2 * ro1 + 1 / (gm2 * R1);       % denomerator
gain = R1 * (Vout ./ Vn2);                                                      % Gain: Vn_out/Vn2

% ===========
% Plots
% ===========
figure;
semilogx(frequencies, 20 * log10(abs(gain)), 'LineWidth', 2); % Gain in dB
grid on;
title('Frequency Response of Cascode Common Gate (CG)', 'FontSize', 20);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Gain (dB)', 'FontSize', 20);
legend('|V_{n,out}/V_{n2}|', 'FontSize', 20);
