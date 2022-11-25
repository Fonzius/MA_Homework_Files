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
title('Input impedance of exponential horn with radiation impedance')

% Calculate input impedance Z1_ori without radiation
Z_end_0=zeros(1,length(f));
for i=length(x):-1:2
    Zin1_ori=ZIN1(a(i-1),a(i),delta_ultimate,Z_end_0,k,rho,c);
    Z_end_0=Zin1_ori;
end

figure();
plot(f,20*log10(abs(Zin1_ori)))
title('Input impedance of exponential horn without radiation impedance')


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
Zin_compound_d=1j.*Zin_exp.*sin(k.*L_cyl)+Z0.*cos(k.*L_cyl);
Zin_compound=Z0.*Zin_compound_n./Zin_compound_d;


[Zin1_ori_p,f1_ori_p]=findpeaks(20*log10(abs(Zin1_ori)),f);
[Zin1_rad_p,f1_rad_p]=findpeaks(20*log10(abs(Zin1_rad)),f);
[Zin1_compound_p,f1_compound_p]=findpeaks(20*log10(abs(Zin_compound)),f);



figure();
plot(f,20*log10(abs(Zin_compound)))
title('Input impedance of compound horn with radiation impedance')

figure();
subplot(2,1,1)
plot(f,20*log10(abs(Zin1_ori)),'k',f,20*log10(abs(Zin1_rad)),'b',...
    f,20*log10(abs(Zin_compound)),'r',f1_ori_p,Zin1_ori_p,'k*',...
    f1_rad_p,Zin1_rad_p,'b*',f1_compound_p,Zin1_compound_p,'r*')
legend('exponential-without radiation','exponential-with radiation','compound')

subplot(2,1,2)
plot(f,angle(Zin1_ori),'k',f,angle(Zin1_rad),'b',f,angle(Zin_compound),'r')
legend('exponential-without radiation','exponential-with radiation','compound')

x_ori=1:length(f1_ori_p);
x_rad=1:length(f1_rad_p);
x_compound=1:length(f1_compound_p);

f_ori_harm=x_ori.*f1_ori_p(1);
f_rad_harm=x_rad.*f1_rad_p(1);
f_compound_harm=x_compound.*f1_compound_p(1);

figure();
subplot(3,1,1)
scatter(x_ori,f1_ori_p,'k')
hold on
scatter(x_ori,f_ori_harm,'r')
legend('Z-max','Z-harmonic')


subplot(3,1,2)
scatter(x_rad,f1_rad_p,'k')
hold on
scatter(x_rad,f_rad_harm,'r')
legend('Z-max','Z-harmonic')

subplot(3,1,3)
scatter(x_compound,f1_compound_p,'k')
hold on
scatter(x_compound,f_compound_harm,'r')
legend('Z-max','Z-harmonic')


