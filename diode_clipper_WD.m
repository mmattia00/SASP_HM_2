clear; close all; clc

%% Simscape File for Ground Truth

load('ssc_output.mat')

%% Sampling Frequency
fs = 96e3;

%% Sampling Period
Ts = 1/fs;

%% Simulation Duration
stop_time = 1;  % [seconds]

%% Input Signal
% Fundamental Frequency
f0 = 440;

% Time Axis
t = 0:Ts:stop_time;

% Signal Amplitude
A = 1.5;
vin = A * sin(2*pi*f0*t);

% Music Signal Test (Uncomment the following lines for test)
% [vin, fs_rec] = audioread('guitar_input.wav');
% G = 5;                                  % Gain factor 
% vin = G * (vin(:, 1) + vin(:, 2))/2;    % Convert the signal to MONO
% vin = resample(vin, fs, fs_rec);        % Resampling the signal from 44.1 kHz to 96 kHz
% t = 0:Ts:Ts*(length(vin)-1);

%% Circuit Parameters
% Resistive Elements
Rin = 1;
R1 = 1e4;
Rout = 1e4;

% Dynamic Elements
C1 = 1e-6;
C2 = 1e-9;

%% Setting of Free Parameters (Adaptation Conditions)
Z5 = Rout;
Z6 = Ts/(2*C2);
Z9 = R1;
Z12 = Ts / (2*C1);
Z11 = Rin;
% Makes ports of adaptors reflection-free
Z4 = (Z5*Z6)/(Z5+Z6);
Z1 = Z4;
Z10 = Z11 + Z12;
Z8 = Z10;
Z7 = Z8 + Z9;
Z2 = Z7;
Z3 = (Z1*Z2)/(Z1 + Z2);

%% Computing Scattering Matrices
Bser = [1, 1, 1];
Qpar = [1, 1, 1];
Zser1 = diag([Z10, Z11, Z12]);
Zser2 = diag([Z7, Z8, Z9]);
Zpar1 = diag([Z1, Z2, Z3]);
Zpar2 = diag([Z4, Z5, Z6]);

Sser1 = eye(3) - 2*Zser1*Bser'*inv(Bser*Zser1*Bser')*Bser;
Sser2 = eye(3) - 2*Zser2*Bser'*inv(Bser*Zser2*Bser')*Bser;
Spar1 = 2*Qpar' * inv(Qpar*inv(Zpar1)*Qpar')*Qpar*inv(Zpar1) - eye(3);
Spar2 = 2*Qpar' * inv(Qpar*inv(Zpar2)*Qpar')*Qpar*inv(Zpar2) - eye(3);



%% Initialization of Waves
b12 = 0;
b6 = 0;
a5 = 0;
a9 = 0;

%% Initialization of Output Signals
vout = zeros(1, length(t));

%% Simulation Algorithm

for n = 1 : length(t)

    % input
    a11 = vin(n);

    % dynamic elements
    a12 = b12;
    a6 = b6;
    
    % FORWARD
    b4 = Spar2(1, :)*[0; a5; a6];
    a1 = b4;
    
    b10 = Sser1(1, :)*[0; a11; a12];
    a8 = b10;
    
    b7 = Sser2(1, :)*[0; a8; a9];
    a2 = b7;

    b3 = Spar1(3, :)*[a1; a2; 0];

    % non-linear component
    a3 = antiparallel_diodes(b3, Z3);

    % BACKWARD
    b1 = Spar1(1, :)*[a1; a2; a3];
    a4 = b1;

    b2 = Spar1(2, :)*[a1; a2; a3];
    a7 = b2;

    b6 = Spar2(3, :)*[a4; a5; a6];
    
    b8 = Sser2(2, :)*[a7; a8; a9];

    a10 = b8;

    b12 = Sser1(3, :)*[a10; a11; a12];
    % useless to compute the 2 reflected waves of resistences 
    % b11 = Sser1(2, :)*[a10; a11; a12];
    % b9 = Sser2(3, :)*[a7; a8; a9];
    b5 = Spar2(2, :)*[a4; a5; a6];
   
    % Read Output
    vout(n) = (a5+b5)/2;   
end

% Uncomment the following line to hear the Diode Clipper
% sound(vout, fs)

%% Output Plots

plot_lim = 5/f0; % Limit the plot to just 5 periods of the output signal

figure
set(gcf, 'Color', 'w');
plot(gt(1, :), gt(2, :), 'r', 'Linewidth', 2);
hold on;
plot(t, vout, 'b--', 'Linewidth', 2);
grid on;
xlim([0, plot_lim]);
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{out}}$ [V]','Fontsize',16,'interpreter','latex');
legend('Simscape','WDF','Fontsize',16,'interpreter','latex');
title('Output Signal','Fontsize',18,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, vout - gt(2, :), 'k', 'Linewidth', 2);
grid on;
xlim([0, plot_lim]);
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{out}}$ [V]','Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',18,'interpreter','latex');



%% Compute Mean Squared Error (MSE)

mse = mean((vout - gt(2, :)).^2);
disp('MSE = ')
disp(mse)

%% Compute spectrum to visualize the harmonics introduced by distortion

frequency_axis = fs/2*linspace(-1,1,fs+1);
figure;
plot(frequency_axis, 20*log(abs(fft(vin))));
hold on;
plot(frequency_axis, 20*log(abs(fft(vout))));
title("Magnitude FFT of input vs output")
xlabel('Frequency (Hz)');
ylabel('magnitude');
legend("input signal", "output signal");

