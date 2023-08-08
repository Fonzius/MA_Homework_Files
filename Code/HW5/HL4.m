clear all;
close all;
clc;
%%
close all;

L = 1.37;
r = 0.6e-2;
c = 343.2816;
f = linspace(50,1200,256);
k = 2*pi*f/c;
rho = 1.2039;

%%
Z_trumpet = importdata("trumpet.txt");
for ii = 1:3
    Z_trumpet(1,:) = [];
end
% figure()
Z_trumpet = str2double(Z_trumpet);
 plot(f,Z_trumpet(:,2)')
figure()
[~, fpeakbell_trumpet]= findpeaks(Z_trumpet(:,2)',f);


%%
Z_bell = importdata("tube with bell.txt");
for ii = 1:3
    Z_bell(1,:) = [];
end
% figure()
Z_bell = str2double(Z_bell);
%  plot(f,Z_bell(:,2)')
figure()
[~, fpeakbell]= findpeaks(Z_bell(:,2)',f);

% Curve fitting
intervalfreq= diff(fpeakbell);
%%%% Method 1
averintervalfreq=sum(intervalfreq)/length(intervalfreq);
multipule1=fpeakbell(1)/averintervalfreq*1.1;
multipule1=0.55;

xpeakbell=linspace(1,length(fpeakbell),length(fpeakbell));
xpeakbell_=linspace(1,length(fpeakbell),length(fpeakbell));
xpeakbell_(1)=[];
xpeakbell_=[multipule1 xpeakbell_];
ystandardbell=xpeakbell_./xpeakbell_(1).*fpeakbell(1);
% ystandardbell2=fpeakbell(1).*xpeakbell;

%%%% Method 2
% first fitting the data without the first one with given slope 1
% pfit_x=xpeakbell;
% pfit_x(end)=[];

pfit_x=0:length(xpeakbell)-1;
fittotal=polyfit(pfit_x,fpeakbell,1);
fitslope=fittotal(1);
fitb=fittotal(2);
yfit=fitslope.*pfit_x+fitb;
% scatter(pfit_x,fpeakbell)
% hold on
% plot(pfit_x,yfit)
% legend('raw','fit')

intervalfreq_fit=diff(yfit);
averintervalfreq2=intervalfreq_fit(1);
multiple2=fpeakbell(1)/averintervalfreq2;

xpeakbell_2=linspace(1,length(fpeakbell),length(fpeakbell));
xpeakbell_2(1)=[];
xpeakbell_2=[multiple2 xpeakbell_2];
ystandardbell2=xpeakbell_2./xpeakbell_2(1).*fpeakbell(1);


freq_fit=fpeakbell(1)./multiple2.*xpeakbell_2(2:end);


scatter(xpeakbell,fpeakbell)
% hold on
% plot(xpeakbell,ystandardbell) % method 1
hold on
plot(xpeakbell,[xpeakbell_2(1) freq_fit]) % method 2
% set(gca,'xtick',xpeakbell_) 
legend({'raw','method1','method2'},"Location","southeast")
xlabel('Mutiples')
ylabel('Frequency(Hz)')

%%

Z_comsol = importdata("tube.txt");
for ii = 1:3
    Z_comsol(1,:) = [];
end
Z_comsol = str2double(Z_comsol);



Z_antonacci = importdata("Tube_Antonacci.txt");
for ii = 1:3
    Z_antonacci(1,:) = [];
end

Z_antonacci = str2double(Z_antonacci);

<<<<<<< HEAD
Zl = (0.25.*(f.*2.*pi).^2.*rho./pi./c+1j.*0.61.*rho.*f*2*pi/pi./r);
Z = ZIN1(r,r,L,Zl,k,rho,c);
figure()
[~, fpeak]= findpeaks(Z_comsol(:,2)',f);
xpeak=linspace(1,length(fpeak),length(fpeak));
xpeakodd=linspace(1,2*length(fpeak)-1,length(fpeak));
ystandard=linspace(1,2*length(fpeak)-1,length(fpeak)).*fpeak(1);
scatter(xpeakodd,fpeak)
hold on
plot(xpeakodd,ystandard)
set(gca,'xtick',xpeakodd) 
legend({'Comsol','Theory'},"Location","southeast")
xlabel('Mutiples')
ylabel('Frequency(Hz)')
figure()

Z0 = (rho*c/(pi*r^2));

Z_L = (0.25.*(f.*2.*pi).^2.*rho./pi./c+1j.*0.61.*rho.*f*2*pi/pi./r);
Z = Z0.* (Zl.*cos(k.*L) + 1j.*Z0.*sin(k.*L))./(1j.*Zl.*sin(k.*L) + Z0.*cos(k.*L));


figure

plot(f,db(Z));

hold on
plot(f,Z_comsol(:,2)', '--');
legend("Matlab Results", "Comsol Simulation");
xlabel("Frequency[Hz]", 'Interpreter','latex');
ylabel("Impedance dB [$\frac{Kg*s}{m^4}$]", 'Interpreter','latex');
grid("minor");
