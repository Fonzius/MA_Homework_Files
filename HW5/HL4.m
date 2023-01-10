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
