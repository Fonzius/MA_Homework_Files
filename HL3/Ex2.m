%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Mirco Pezzoli                                                           %
% 2022   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
%%
addpath('Functions')
simulink_folder = './';       % Simulink projects folder
addpath(simulink_folder);

%% Setup
fs = 44100;                        % Sampling frequency
signalLen = 3;                     % Signal length
t = 0 : 1/fs : signalLen-1/fs;     % Time axis

fileName = 'complete_guitar.wav';  % Audio file path

%% Simulation
% run the simulink simulation using the command sim (see doc sim).
sim("Ex2.slx");
%%
% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1

% We used a fixed time interval of 1/44100 in the simulink file so
% resampling is not necessary since all the samples are equally spaced in
% time


% Plot the resampled signal in time
figure()
plot(I1.Time, I1.Data, LineWidth=2);
ylim([-1e-3 1e-3]);
xlabel("Time [s]");
ylabel("Current [A]");
%title("Guitar Signal");

figure()
plot(noString.Time, noString.Data, LineWidth=2);
xlabel("Time [s]");
ylabel("Current [A]");
%ylim([-1e-3 1e-3]);
%title("Guitar Signal Nostring");


% Normalize the signal
sound = (I1.Data./max(abs(I1.Data)));

%% Plot and play
% Plot the signal frequency content as magnitude and phase
sound_fft = myFFT(sound,fs);
f = 1:1:length(sound_fft);

figure(1)
plotFFT_linearFreqScale(abs(sound_fft), angle(sound_fft), f, f(2)-f(1), fs, fs/2,1);


sound_nostring = (noString.Data./max(abs(noString.Data)));
sound_fft_nostring = myFFT(sound_nostring,fs);

figure(2)
plotFFT_linearFreqScale(abs(sound_fft_nostring), angle(sound_fft_nostring), f, f(2)-f(1), fs, fs/2,2);


% Play the sound
soundsc(sound, fs);

% Save on disk
audiowrite('10669941_Bernasconi_10876787_Luan_Guitar.wav', sound, fs);
audiowrite('10669941_Bernasconi_10876787_Luan_Guitar_NoString.wav', sound_nostring, fs);

disp('Save file on disk...')
%%
winwid=0.04; % 窗宽度
at_thd=5e-3; %起始时刻阈值
pw_ref=1e-6;
dt =1/fs;
Nw=round(winwid/dt);
Nov=round(Nw/5);



[sample, FS] = audioread("GuitarSample.wav");

sample_fft = myFFT(sample(:,1), fs);

figure(1)
spectrogram(sound', Nw, Nov, [], fs);
xlim([0,8]);
ylim([0 3]);
title("Spectrogram of Guitar simulation with string simulation");
figure(2)
spectrogram(sound_nostring, Nw, Nov, [], fs);
xlim([0,8]);
ylim([0 3]);
title("Spectrogram of Guitar simulation without string simulation");
figure(3)
spectrogram(sample(:,1), Nw, Nov, [], fs);
xlim([0,8]);
ylim([0 3]);
title("Spectrogram of a real Guitar");
