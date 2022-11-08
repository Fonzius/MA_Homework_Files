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

set_param('HL2/R', 'R', num2str(R));
set_param('HL2/C', 'C', num2str(C));
set_param('HL2/L', 'L', num2str(L));

alpha = R/(2*L);
%% Ex 1.b
f = 0:Fs/length(out.impulse.Data):Fs-(1/length(out.impulse.Data));
H = (abs(fft(out.impulse_response.Data)./(fft(out.impulse.Data))));
 
figure()
plot(f,H);
xlim([0,Fs/2]);

f_0 = c/(2*pi) * sqrt(S/(V*l));

%% Ex 2.a
set_param('HL2_ex2/R', 'R', num2str(R));
set_param('HL2_ex2/C', 'C', num2str(C));
set_param('HL2_ex2/L', 'L', num2str(L));

H_branch = (abs(fft(out.output.Data)./(fft(out.input.Data))));
 
figure()
plot(f,H_branch);
xlim([0,Fs/2]);