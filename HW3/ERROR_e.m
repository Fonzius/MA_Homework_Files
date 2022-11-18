function [e1, delta]= ERROR_e (n_L)
% HW3
% clc
% clear
% clf

L=0.35; % [m]
a0=0.008; % [m]
m=4.2; % [m^-1]
rho=1.2; % [kg/m^3]  (20[c])
c=343; % [m/3] (20[c])

fmax = 2000;          % maximum evaluation frequency (Hz)
N = fmax;             % number of frequencies for evaluation (even)
finc = fmax / (N-1);
f = 1:finc:fmax; % eps  2.2204e-5
omega=2*pi*f;
k=omega/c;
%//////WITHOUT WALL LOSSES///////

% Set up geometry

% n_L=100;
delta=L/n_L; %length interval [m]
x=0:delta:L;
x=x';
a=a0.*exp(m.*x);

% Zin=zeros(length(x),length(f));

% Compute Zin - conical segments
Z_end=zeros(1,length(f));
Z_end_0=Z_end;
for i=length(x):-1:2
    Zin1=ZIN1(a(i-1),a(i),delta,Z_end_0,k,rho,c);
    Z_end_0=Zin1;
end

% Compute Zin - exponential horn
b=sqrt(k.^2-m^2);
theta=atan(m./b);
Zin2_n=Z_end.*cos(b.*L+theta)+1j.*(rho.*c./(pi.*a(end)^2)).*sin(b.*L);
Zin2_d=1j.*Z_end.*sin(b.*L)+(rho.*c./(pi.*a(end).^2)).*cos(b.*L-theta);
Zin2=rho.*c./(pi*a(1).^2).*Zin2_n./Zin2_d;

figure();
% subplot(3,1,1)
plot(f,20*log10(abs(Zin1)),'k',f,20*log10(abs(Zin2)),'r--','LineWidth',1)
title('abs')
legend('Z1-conical','Z2-exponential')

% subplot(3,1,2)
% plot(f,20*log10(real(Zin1)),'k',f,20*log10(real(Zin2)),'r--','LineWidth',1)
% title('real')
% legend('Z1','Z2')
% 
% subplot(3,1,3)
% plot(f,20*log10(imag(Zin1)),'k',f,20*log10(imag(Zin2)),'r--','LineWidth',1)
% title('imag')
% legend('Z1','Z2')

%% a)

e_eff=1./(omega(end)-omega(1));
Zin_delta=(abs(Zin1-Zin2)).^2;
e1=e_eff.*trapz(Zin_delta);

%% b)

figure();
findpeaks(abs(Zin1),f);
[Zin1_p,f1_p]=findpeaks(abs(Zin1),f);

rr=11;