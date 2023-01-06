clear all;
close all;
clc;
%%
close all;

L = 1.33;
r = 0.6e-2;
c = 343;
f = linspace(50,1200,256);
k = 2.*pi.*f/c;
rho = 1.225;


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

Zl = (0.25.*(f.*2.*pi).^2.*rho./pi./c+1j.*0.61.*rho.*f*2*pi/pi./r);
Z = ZIN1(r,r,L,Zl,k,rho,c);

figure
plot(f,db(Z));
hold on
plot(f,Z_comsol(:,2)');
hold on
plot(f,Z_antonacci(:,2)','--');
legend("Matlab", "Ours", "Antonacci");
