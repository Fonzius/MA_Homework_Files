clear all;
clc;
%% Data
h = 0.25;
w = 0.25;
l = 0.18;

ln = 0.06;
rn = 0.025;

c = 343;

%% Point a
f_0 = c/(2*pi) * sqrt(rn^2*pi/(h*w*l*ln));

%% Point b
ln_1 = ln + 16*rn/(3*pi);
f_0_v = c/(2*pi) * sqrt(rn^2*pi/(h*w*l*ln_1));

%% Point c
rho_air = 1.2754;
omega_0_v = 2*pi*f_0_v;
m = rho_air*ln_1*pi*(rn^2);
R = 2*m*omega_0_v; %[Kg/s]

%% Point d
R_1 = 5*10^(-4);
k = omega_0_v^2 * m;
freq = 0:1:200;
omega = 2*pi*freq;

Z = R_1 + 1i*(omega.*m - k./omega);

figure();
subplot(2,1,1);
semilogy(freq,abs(Z));
title("Absolute value of Impedance");
ylabel("Module [Kg/s]");
xlabel("f [hz]");
subplot(2,1,2);
semilogy(freq,angle(Z));
title("Phase of Impedance")
ylabel("Phase [rad/s]");
xlabel("f [hz]");

%% Point e
alpha = R_1/(2*m);

Q = omega_0_v / (2*alpha);
tau = 1/alpha;

%% Point f
R_2 = 0:0.01:0.5;
alpha_2 = R_2/(2*m);
Q_2 = omega_0_v ./ (2.*alpha_2);

figure();
plot(R_2,Q_2);


