% HW3
clc
clear

L=0.35; % [m]
a0=0.008; % [m]
m=4.2; % [m^-1]
rho=1.2; % [kg/m^3]  (20[c])
c=343; % [m/3] (20[c])
f=0:1:2000;
omega=2*pi*f;
k=omega/c;

% Set up geometry

n_L=100;
delta=L/n_L; %length interval [m]
x=0:delta:L;
x=x';
a=a0.*exp(m.*x);

% Zin=zeros(length(x),length(f));

% Compute Zin - conical segments
Z_end=zeros(1,length(f));
for i=length(x):-1:2
    Zin1=ZIN(i-1,i,delta,Z_end,k,rho,c);
    Z_end=Zin1;
end



% Compute Zin - exponential horn
b=sqrt(k.^2-m^2);
theta=atan(m./b);
Zin2_n=Z_end.*cos(b.*L+theta)+1j.*(rho.*c./(pi.*a(end)^2)).*sin(b.*L);
Zin2_d=1j.*Z_end.*sin(b.*L)+(rho.*c./((pi.*a(end)^2).*cos(b.*L-theta)));
Zin2=rho.*c./(pi*a(1)^2).*Zin2_n./Zin2_d;


% plot(f,abs(Zin1),'b',f,abs(Zin2),'r')
%plot(f,abs(Zin2))




