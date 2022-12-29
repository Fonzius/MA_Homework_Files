clc;
clear all;
close all;

%% Point b
L = 1;
H = 1.4;

%PlateZ: old results, very high impedance, PlateZ2: new results, good
%impedance, bad results

% Z = importdata("PlateZ.txt");
% 
% for ii=1:3 % Removing empty rows
%     Z(1,:) = [];
% end

Z = importdata("PlateZ2.txt");

for ii=1:4 % Removing empty rows
    Z(1,:) = [];
end

% legend of column indices and their meaning
% 1 = x
% 2 = y
% 3 = f
% 4 = abs(Z)
% 5 = Z

% Sort the whole array of Z by frequency
Z_sorted = sortrows(Z,3);

%Sort for each frequency by db(abs(Z)) in a separate array
Z_f1 = str2double(sortrows(Z_sorted(1:48,:),4));
Z_f2 = str2double(sortrows(Z_sorted(49:96,:),4));
Z_f3 = str2double(sortrows(Z_sorted(97:144,:),4));
Z_f4 = str2double(sortrows(Z_sorted(145:192,:),4));
Z_f5 = str2double(sortrows(Z_sorted(193:240,:),4));


% Plotting the first 6-7 points with the lowest abs(Z) to do the design by
% hand by choosing the individual points where to put the strings

figure
plot(Z_f1(1:10,1),Z_f1(1:10,2), '-o', Color='r');
hold on
plot(Z_f2(1:10,1),Z_f2(1:10,2), '-o', Color='g');
hold on
plot(Z_f3(1:10,1),Z_f3(1:10,2), '-o', Color='b');
hold on
plot(Z_f4(1:10,1),Z_f4(1:10,2), '-o', Color='k');
hold on
plot(Z_f5(1:10,1),Z_f5(1:10,2), '-o', Color='magenta');
hold on

plot(L-Z_f1(1:10,1),Z_f1(1:10,2), '-o', Color='r');
hold on
plot(L-Z_f2(1:10,1),Z_f2(1:10,2), '-o', Color='g');
hold on
plot(L-Z_f3(1:10,1),Z_f3(1:10,2), '-o', Color='b');
hold on
plot(L-Z_f4(1:10,1),Z_f4(1:10,2), '-o', Color='k');
hold on
plot(L-Z_f5(1:10,1),Z_f5(1:10,2), '-o', Color='magenta');


%Chosen points for the bridge. Values are for PlateZ so the bridge for
%PlateZ2 will look weird

%For PlateZ
% Z_final = [Z_f1(2,1), H-Z_f1(2,2), Z_f1(2,3:5); ...
%            Z_f2(1,1), H-Z_f2(1,2), Z_f2(1,3:5); ...
%            Z_f3(6,1), H-Z_f3(6,2), Z_f3(6,3:5); ...
%            L-Z_f4(2,1), Z_f4(2,2), Z_f4(2,3:5); ...
%            L-Z_f5(6,1), Z_f5(6,1), Z_f5(6,3:5); ...
%            ];

% For PlateZ2
Z_final = [Z_f1(6,1), H-Z_f1(6,2), Z_f1(6,3:5); ...
           Z_f2(1,1), H-Z_f2(1,2), Z_f2(1,3:5); ...
           Z_f3(1,1), H-Z_f3(1,2), Z_f3(1,3:5); ...
           L-Z_f4(2,1), Z_f4(2,2), Z_f4(2,3:5); ...
           L-Z_f5(4,1), Z_f5(4,1), Z_f5(4,3:5); ...
           ];
figure
plot(Z_final(:,1), Z_final(:,2));
axis([0,L, 0,H]);
title("Bridge Design");

%% Point c
rho = 10.8e-3; %[kg/m] piano wire linear density, from the HW text

f = [349.23, 440, 523.25, 659.25, 783.99];
omega = f.*2.*pi;

c = 2.*Z_final(:,2).*f'; %lambda * f
T = rho.*c.^2;

Z0 = T./c;
Y = 1./(Z_final(:,5));

% Acoustics of musical instruments, page 275, eq 6.34 (also on Antonacci's slides)
% This is the source of all the problems probably. It's way too small for
% PlateZ and since real(X) is what determines the damping the string will
% basically not dampen. For PlateZ2 for some reason its values will break
% everything that comes after this


X = (Y .* 1j .* Z0)./pi;

a = zeros(length(f), 300);
b = zeros(length(f), 300);
epsilon = zeros(length(f), 300);

% Detune of the string. Built so that omega*epsilon always goes from -3 to 3 
% as in the graphs
for ii = 1:length(f)
    epsilon(ii,:) = linspace(-3/(omega(ii)), 3/(omega(ii)), 300);
end

% Acoustics of musical instruments, page 279, eq 6.50 (also on Antonacci's slides)
for ii = 1:length(f)

    a(ii,:) = X(ii) + epsilon(ii,:) + sqrt(X(ii).^2 + epsilon(ii,:).^2);
    b(ii,:) = X(ii) + epsilon(ii,:) - sqrt(X(ii).^2 + epsilon(ii,:).^2);
    
end

%Plot of eigenfrequency shift wrt string detune. Should be like graphs at
%page 279. For values in PlateZ the graph looks convincing, for PlateZ2 not
%as much: either a or b is way too big or way too small.

figure
for ii = 1:length(f)
    subplot(length(f),1,ii);
    plot(omega(ii).*epsilon(ii,:), omega(ii).*real(a(ii,:)));
    hold on
    plot(omega(ii).*epsilon(ii,:), omega(ii).*real(b(ii,:)));
    title(ii + "th freq " + f(ii));
end

%% Point d

t = linspace(0,10,1000); %time axis

x = zeros(length(f), 1000); %x axis, goes from 0 to the length (Z_final(ii,2)) of the string
for ii = 1:length(f)
    x(ii,:) = linspace(0,Z_final(ii,2),1000);
end

%Displacement computed for an arbitrary epsilon. The higher/lower epsilon
%should give more minima, epsilon values around the middle (150th index)
%should give less minima. For PlateZ if you increase by multiple orders of
%magnitude X the reults are very convincing, for PlateZ2 it seems that the
%damping is way too strong or something

F0 = 1; % F(0), initial condition for force on the strings
mu = sqrt(epsilon(:,10).^2 + (X.^2)); 

% Acoustics of musical instruments, page 280
Vb = ((2.*pi.*F0.*X)./(mu.*Z0)) .* exp(1j.*(epsilon(:,10) + X + 1) .* (omega'.*t)) .* (mu.*cos(mu.*omega'.*t) + 1j.*X.*sin(mu.*omega'.*t));

figure
plot(t, db(abs(Vb(5,:)./Vb(5,1))));
title("Bridge displacement");

