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

set_param('HL2_ex1', 'PreLoadFcn', num2str(Fs));
set_param('HL2_ex1/Resonator/R', 'R', num2str(R));
set_param('HL2_ex1/Resonator/C', 'C', num2str(C));
set_param('HL2_ex1/Resonator/L', 'L', num2str(L));

out = sim('HL2_ex1', dur);

%% Ex 1.b
f = 0:Fs/length(out.impulse.Data):Fs-(1/length(out.impulse.Data));
H = (abs(fft(out.impulse_response.Data)./(fft(out.impulse.Data))));
 
figure()
plot(f,H);
xlim([0,Fs/2]);

f_0 = c/(2*pi) * sqrt(S/(V*l));

%% Ex 2.a

subsystems = find_system('HL2_ex2a','SearchDepth',2,'BlockType','SubSystem');

for ii = 1:length(subsystems)
    if contains(subsystems(ii), 'Resonator')
        set_param(char(append(subsystems(ii),'/R')), 'R', num2str(R));
        set_param(char(append(subsystems(ii),'/C')), 'C', num2str(C));
        set_param(char(append(subsystems(ii),'/L')), 'L', num2str(L));
    end
end
set_param('HL2_ex2a', 'PreLoadFcn', num2str(Fs));
%%
H_2a = (abs(fft(out.output.Data)./(fft(out.input.Data))));

figure()
plot(f,H_2a);
xlim([0,Fs/2]);
title("2 branch");
%% Ex 2.b
subsystems = find_system('HL2_ex2b','SearchDepth',2,'BlockType','SubSystem');

for ii = 1:length(subsystems)
    if contains(subsystems(ii), 'Resonator')
        set_param(char(append(subsystems(ii),'/R')), 'R', num2str(R));
        set_param(char(append(subsystems(ii),'/C')), 'C', num2str(C));
        set_param(char(append(subsystems(ii),'/L')), 'L', num2str(L));
    end
end
set_param('HL2_ex2b', 'PreLoadFcn', num2str(Fs));
%%
H_2b = (abs(fft(out.output.Data)./(fft(out.input.Data))));

figure()
plot(f,H_2b);
xlim([0,Fs/2]);
title("3 branch");