clc;
clear all;
close all;

alpha = 0.75 * pi/180; %[rad]
L = 0.45; %[m]
f_e4 = 329.63; %[hz]
c = 343; %[m/s]
rho = 1.225;

L_c = 0.04; %[m] Mouth end correction for low freqs

%% Question 1

r0 = linspace(L*tan(alpha),1, 10000);
r1 = r0 - L*tan(alpha);

k_e4 = 2*pi*f_e4/c;

Z_mouth = 1j* (L_c*rho./(pi.*r0.^2)) .* 2.*pi.*f_e4;



a0 = r0(min_index_1);
a1 = r1(min_index_1);

x1 = a1/tan(alpha);

