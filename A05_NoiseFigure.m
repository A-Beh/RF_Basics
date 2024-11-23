% Ali Behfarnia
% Edited 11/2024
% Goal: Simple Example of the Noise Figure (NF)
% Three NF elements: 
% One for noise antenna, One for BPF (to supress out of band), and One for LNA

clc; clear all; close all;

% ===========
% Noise Figures in dB
% ===========
NF_antenna_dB = 3; % Antenna
NF_BPF_dB = 2;     % Band-pass filter
NF_LNA_dB = 1;     % Low-noise amplifier

% ===========
% Available Power Gain in dB
% ===========
G_BPF_dB = 0;       % BPF
G_LNA_dB = 20;     % LNA

% ===========
% dB to Linear 
% ===========
NF_antenna = 10^(NF_antenna_dB / 10);
NF_BPF = 10^(NF_BPF_dB / 10);
NF_LNA = 10^(NF_LNA_dB / 10);

G_BPF = 10^(G_BPF_dB / 10);
G_LNA = 10^(G_LNA_dB / 10);

% ===========
% Friis Equation: Total NF
% ===========
NF_total = NF_antenna + (NF_BPF - 1) / G_BPF + (NF_LNA - 1) / (G_BPF * G_LNA);
NF_total_dB = 10 * log10(NF_total);

% ===========
% Results
% ===========
disp(['Total Noise Figure (NF_total) in linear scale: ', num2str(NF_total)]);
disp(['Total Noise Figure (NF_total) in dB: ', num2str(NF_total_dB), ' dB']);
