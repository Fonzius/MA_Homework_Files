clc;
clear all;
close all;

%% Ex 1.a
V = 0.1; %[m^3]
l = 0.1; %[m]
S = 100; %[m^2]
c = 343; %[m/s]
rho = 1.2; %[Kg/m^3]

Fs = 50000;

R = rho*c/S;
L = rho*l/S;
C = V/(rho*c^2);

%% Setup all the system R,L,C
subsystems = find_system('HL2','SearchDepth',2,'BlockType','SubSystem');

for ii = 1:length(subsystems)
    if contains(subsystems(ii), 'Resonator')
        set_param(char(append(subsystems(ii),'/R')), 'R', num2str(R));
        set_param(char(append(subsystems(ii),'/C')), 'C', num2str(C));
        set_param(char(append(subsystems(ii),'/L')), 'L', num2str(L));
    end
end

set_param('HL2', 'PreLoadFcn', num2str(Fs));

% Run the simulation before going on with the next parts of the code

%% Ex 1.b
f = 0:Fs/length(out.in1x1.Data):Fs-(1/length(out.in1x1.Data));
H = (abs(fft(out.out1x1.Data)./(fft(out.in1x1.Data))));

f_0 = c/(2*pi) * sqrt(S/(V*l));

maxh = 0;
maxi = 0;
for ii = 1:length(H)
    if H(ii)>maxh
        maxh = H(ii);
        maxi = ii;
    end
end

figure()
plot(f,H,'LineWidth',2);
xlim([0,Fs/2]);
xline(f_0, 'b', 'LineWidth',2);
xline(f(maxi), '--r', 'LineWidth',2);
ylabel("Admittance [S]");
xlabel("f [Hz]");
title("Frequency Response of the Electric-Equivalent of an Helmholtz Resonator")

error = (f(maxi)-f_0)/100;

%% Ex 2.a

H_2 = (abs(fft(out.out2x2.Data)./(fft(out.in2x2.Data))));

figure()
plot(f,H_2);
xlim([0,Fs/2]);
title("2x2");

%% Ex 2.b

H_2b = (abs(fft(out.out3x3.Data)./(fft(out.in3x3.Data))));

figure()
plot(f,H_2b);
xlim([0,Fs/2]);
title("3 branch");