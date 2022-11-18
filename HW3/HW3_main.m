% HW3_main
clc
clear


%% a)
n_L=1:10;%length interval [m] 1:100
e1=zeros(length(n_L),1);
e2=zeros(length(n_L),1);
delta=zeros(length(n_L),1);

for i=1:length(n_L)
    [e1(i,1),e2(i,1),delta(i,1)]=ERROR_e(n_L(i));
end
figure();
plot(delta,e1)
title('e1')

%% b)
figure();
% plot(delta,e2)
scatter(delta,e2)
title('e2')

% delta_ultimate=delta(6); %0.0583

%% c)
delta_ultimate=0.0583;
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

x=0:delta_ultimate:L;
x=x';
a=a0.*exp(m.*x);
Z_L0=0.25.*omega.^2.*rho./pi./c+1j.*0.61.*rho.*omega./pi./a(end);
S_p=pi*a(end)^2;
theta= atan(delta_ultimate/(a(end)-a(end-1))); % flaring angle of the last conical section
S_s=2*S_p/(1+cos(theta));
Z_L=Z_L0.*S_p./S_s;

% Calculate input impedance Z1_rad with radiation
Z_end_0=Z_L;
for i=length(x):-1:2
    Zin1_rad=ZIN1(a(i-1),a(i),delta_ultimate,Z_end_0,k,rho,c);
    Z_end_0=Zin1_rad;
end

figure();
plot(f,20*log10(abs(Zin1_rad)))
title('Input impedance with radiation')
% Calculate input impedance Z1_ori without radiation
Z_end_0=zeros(1,length(f));
for i=length(x):-1:2
    Zin1_ori=ZIN1(a(i-1),a(i),delta_ultimate,Z_end_0,k,rho,c);
    Z_end_0=Zin1_ori;
end

e_eff=1./(omega(end)-omega(1));
Zin_delta=(abs(Zin1_ori-Zin1_rad)).^2;
e1=e_eff.*trapz(Zin_delta);

%% d)
Zin_exp=Zin1_rad;
a_cyl=a(1);
L_cyl=0.85;
S_cyl=pi*a_cyl^2;
Z0=rho*c/S_cyl;

Zin_compound_n=Zin_exp.*cos(k.*L_cyl)+1j.*Z0.*sin(k.*L_cyl);
Zin_compound_d=1j.*Zin_exp.*sin(k.*L)+Z0.*cos(k.*L_cyl);
Zin_compound=Z0.*Zin_compound_n./Zin_compound_d;

figure();
plot(f,20*log10(abs(Zin_compound)))
title('Input impedance of compound horn')


