% HW3_main
clc
clear
close all


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
xlabel('\delta [m]')
ylabel('e1')


%% b)
figure();
% plot(delta,e2)
sz = 25;
ccc = linspace(1,10,length(delta));
scatter(delta,e2,sz,ccc,'filled')
xlabel('\delta [m]')
ylabel('e2')

% delta_ultimate=delta(6); %0.0583

%% c)
delta_ultimate=0.0583;
L=0.35; % [m]
a0=0.008; % [m]
m=4.2; % [m^-1]
rho=1.2; % [kg/m^3]  (20[c])
c=343; % [m/3] (20[c])

fmax = 20000; %2000         % maximum evaluation frequency (Hz)
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

% figure();
% plot(omega,20*log10(abs(Zin1_rad)),LineWidth=1.5)
% xlabel('\omega [rad/s]')
% ylabel('Z_1 [dB]')
% title('Input impedance of exponential horn with radiation impedance')

% Calculate input impedance Z1_ori without radiation
Z_end_0=zeros(1,length(f));
for i=length(x):-1:2
    Zin1_ori=ZIN1(a(i-1),a(i),delta_ultimate,Z_end_0,k,rho,c);
    Z_end_0=Zin1_ori;
end

figure();
subplot(2,1,1)
plot(f,20*log10(abs(Zin1_rad)),'b',LineWidth=1.5)
hold on

plot(f,20*log10(abs(Zin1_ori)),'r',LineWidth=1.5)
xlabel('f [Hz]')
ylabel('Z_1 Amplitude [dB]')
ylim([50 220])
legend('With radiation','Without radiation')
% title('Input impedance of exponential horn without radiation impedance')

subplot(2,1,2)
plot(f,angle(Zin1_rad),'b',LineWidth=1.5)
hold on

plot(f,angle(Zin1_ori),'r',LineWidth=1.5)
xlabel('f [Hz]')
ylabel('Z_1 Phase')
ylim([-2 3])
legend('With radiation','Without radiation')


e_eff=1./(omega(end)-omega(1));
Zin_delta=(abs(Zin1_ori-Zin1_rad)).^2;
e1=e_eff.*trapz(Zin_delta);

%% d)
Zin_exp=Zin1_rad;
a_cyl=a(1);
L_cyl=0.85; %%%%5
S_cyl=pi*a_cyl^2;
Z0=rho*c/S_cyl;

Zin_compound_n=Zin_exp.*cos(k.*L_cyl)+1j.*Z0.*sin(k.*L_cyl);
Zin_compound_d=1j.*Zin_exp.*sin(k.*L_cyl)+Z0.*cos(k.*L_cyl);
Zin_compound=Z0.*Zin_compound_n./Zin_compound_d;


%%%%%%test cylinder
% Zin_cyl_n=0.*cos(k.*L_cyl)+1j.*Z0.*sin(k.*L_cyl);
% Zin_cyl_d=1j.*0.*sin(k.*L_cyl)+Z0.*cos(k.*L_cyl);
% Zin_cyl=Z0.*Zin_cyl_n./Zin_cyl_d;
% Zin_cyl=1j*Z0.*tan(k.*L_cyl);
% 
% plot(f,20*log10(abs(Zin_cyl)));

%%%%%%%%%
[Zin1_ori_p,f1_ori_p]=findpeaks(20*log10(abs(Zin1_ori)),f);
[Zin1_rad_p,f1_rad_p]=findpeaks(20*log10(abs(Zin1_rad)),f);
[Zin1_compound_p,f1_compound_p]=findpeaks(20*log10(abs(Zin_compound)),f);



figure();
subplot(2,1,1)
plot(f,20*log10(abs(Zin_compound)),LineWidth=1.5)
xlabel('f[Hz]')
ylabel('Amplitude[dB]')
% title('Input impedance of compound horn with radiation impedance')
subplot(2,1,2)
plot(f,angle(Zin_compound),LineWidth=1.5)
xlabel('f[Hz]')
ylabel('Phase')



figure();
subplot(2,1,1)
plot(f,20*log10(abs(Zin1_ori)),'k',f,20*log10(abs(Zin1_rad)),'b',...
    f,20*log10(abs(Zin_compound)),'r',f1_ori_p,Zin1_ori_p,'k*',...
    f1_rad_p,Zin1_rad_p,'b*',f1_compound_p,Zin1_compound_p,'r*',LineWidth=1)
legend('exponential-without radiation','exponential-with radiation','compound')
ylim([50 290])
xlabel('f[Hz]')
ylabel('Amplitude[dB]')

subplot(2,1,2)
plot(f,angle(Zin1_ori),'k',f,angle(Zin1_rad),'b',f,angle(Zin_compound),'r',LineWidth=1)
legend('exponential-without radiation','exponential-with radiation','compound')
xlabel('f[Hz]')
ylabel('Phase')
ylim([-2 3.5])
%% e)
x_ori=1:length(f1_ori_p);
x_rad=1:length(f1_rad_p);
x_compound=1:length(f1_compound_p);

f_ori_harm=x_ori.*f1_ori_p(1);
f_rad_harm=x_rad.*f1_rad_p(1);
f_compound_harm=x_compound.*f1_compound_p(1);

peaknumbers=41;
xx=1:peaknumbers;
yy=1:peaknumbers;
yy2=1:2:2*peaknumbers;


figure();
subplot(2,1,1)
plot(xx,yy,xx,yy2,LineWidth=1.5)
hold on
plot(x_ori,f1_ori_p./f1_ori_p(1),LineWidth=1.5)
hold on
plot(x_rad,f1_rad_p./f1_rad_p(1),LineWidth=1.5)
hold on
% plot(x_rad,f_rad_harm./f_rad_harm(1),'--',LineWidth=1.5)
hold on
plot(x_compound,f1_compound_p./f1_compound_p(1),LineWidth=1.5)
hold on
% plot(x_compound,f_compound_harm./f_compound_harm(1),':',LineWidth=1.5)
legend('Harmonic','Harmonic-odd','Exponential-without radiation',...
    'Exponential-with radiation','Compound-with radiation',Location='northwest')
xlim([1 peaknumbers])
xlabel('Peaknumber')
ylabel('Mulipules')
hold on

subplot(2,1,2)
plot(xx,yy,'--*',xx,yy2,'--*',LineWidth=1.5)
hold on
plot(x_ori,f1_ori_p./f1_ori_p(1),'--*',LineWidth=1.5)
hold on
plot(x_rad,f1_rad_p./f1_rad_p(1),'--*',LineWidth=1.5)
hold on
% plot(x_rad,f_rad_harm./f_rad_harm(1),'--',LineWidth=1.5)
hold on
plot(x_compound,f1_compound_p./f1_compound_p(1),'--*',LineWidth=1.5)
hold on
% plot(x_compound,f_compound_harm./f_compound_harm(1),':',LineWidth=1.5)
legend('Harmonic','Harmonic-odd','Exponential-without radiation',...
    'Exponential-with radiation','Compound-with radiation',Location='northwest')
xlim([1 5])
xlabel('Peaknumber')
ylabel('Mulipules')