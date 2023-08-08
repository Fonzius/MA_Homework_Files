function [x_compound, f_compound]= ChangeLengthhorn (L,peaknumbers)

delta_ultimate=0.0583/5;%0.0583
% L=0.35; % [m]
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

% figure();
% subplot(2,1,1)
% plot(f,20*log10(abs(Zin1_rad)),'b',LineWidth=1.5)
% hold on
% 
% plot(f,20*log10(abs(Zin1_ori)),'r',LineWidth=1.5)
% xlabel('f [Hz]')
% ylabel('Z_1 Amplitude [dB]')
% ylim([50 220])
% legend('With radiation','Without radiation')
% % title('Input impedance of exponential horn without radiation impedance')
% 
% subplot(2,1,2)
% plot(f,angle(Zin1_rad),'b',LineWidth=1.5)
% hold on
% 
% plot(f,angle(Zin1_ori),'r',LineWidth=1.5)
% xlabel('f [Hz]')
% ylabel('Z_1 Phase')
% ylim([-2 3])
% legend('With radiation','Without radiation')
% 

e_eff=1./(omega(end)-omega(1));
Zin_delta=(abs(Zin1_ori-Zin1_rad)).^2;
e1=e_eff.*trapz(Zin_delta);

Zin_exp=Zin1_rad;
a_cyl=a(1);
L_cyl=0.85;  %0.85
S_cyl=pi*a_cyl^2;
Z0=rho*c/S_cyl;

Zin_compound_n=Zin_exp.*cos(k.*L_cyl)+1j.*Z0.*sin(k.*L_cyl);
Zin_compound_d=1j.*Zin_exp.*sin(k.*L_cyl)+Z0.*cos(k.*L_cyl);
Zin_compound=Z0.*Zin_compound_n./Zin_compound_d;


% [Zin1_ori_p,f1_ori_p]=findpeaks(20*log10(abs(Zin1_ori)),f);
% [Zin1_rad_p,f1_rad_p]=findpeaks(20*log10(abs(Zin1_rad)),f);
[Zin1_compound_p,f1_compound_p]=findpeaks(20*log10(abs(Zin_compound)),f);



x_compound=1:length(f1_compound_p);



f_compound=f1_compound_p./f1_compound_p(1);
if length(x_compound )< peaknumbers
    x_zeros=zeros(peaknumbers-length(x_compound));
    x_compound=[x_compound x_zeros];
    f_compound=[F_compound x_zeros];
end

x_compound=x_compound(1:peaknumbers);
f_compound=f_compound(1:peaknumbers);

