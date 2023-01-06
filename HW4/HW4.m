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

Z = importdata("PlateZ22.txt");

for ii=1:4 % Removing empty rows
    Z(1,:) = [];
end

% legend of column indices and their meaning
% 1 = x
% 2 = y
% 3 = f
% 4 = 20*log10(abs(Z))
% 5 = Z

% Sort the whole array of Z by frequency
Z_sorted = sortrows(Z,3);

%Sort for each frequency by db(abs(Z)) in a separate array
Z_f1 = str2double(sortrows(Z_sorted(1:35,:),4));
Z_f2 = str2double(sortrows(Z_sorted(36:70,:),4));
Z_f3 = str2double(sortrows(Z_sorted(71:105,:),4));
Z_f4 = str2double(sortrows(Z_sorted(106:140,:),4));
Z_f5 = str2double(sortrows(Z_sorted(141:175,:),4));


% Plotting the first 6-7 points with the lowest abs(Z) to do the design by
% hand by choosing the individual points where to put the strings

figure
subplot(3,2,1)
plot(Z_f1(1:10,1),Z_f1(1:10,2), '-o', Color='r');
hold on
plot(L-Z_f1(1:10,1),Z_f1(1:10,2), '-o', Color='r');

subplot(3,2,2);
plot(Z_f2(1:10,1),Z_f2(1:10,2), '-o', Color='g');
hold on
plot(L-Z_f2(1:10,1),Z_f2(1:10,2), '-o', Color='g');

subplot(3,2,3);
plot(Z_f3(1:10,1),Z_f3(1:10,2), '-o', Color='b');
hold on
plot(L-Z_f3(1:10,1),Z_f3(1:10,2), '-o', Color='b');

subplot(3,2,4);
plot(Z_f4(1:10,1),Z_f4(1:10,2), '-o', Color='k');
hold on
plot(L-Z_f4(1:10,1),Z_f4(1:10,2), '-o', Color='k');

subplot(3,2,5);
plot(Z_f5(1:10,1),Z_f5(1:10,2), '-o', Color='magenta');
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
Z_final = [Z_f1(2,1), Z_f1(2,2), Z_f1(2,3:5); ...
           Z_f2(1,1), Z_f2(1,2), Z_f2(1,3:5); ...
           Z_f3(5,1), Z_f3(5,2), Z_f3(5,3:5); ...
           L-Z_f4(3,1), Z_f4(3,2), Z_f4(3,3:5); ...
           L-Z_f5(1,1), Z_f5(1,2), Z_f5(1,3:5); ...
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
%T = rho.*c.^2;
T = [1000;1000;1000;1000;1000];

%c = sqrt(T/rho);

%Z0 = T./c;
Z0 = sqrt(T.*rho);

Y = 1./(Z_final(:,5));

% Acoustics of musical instruments, page 275, eq 6.34 (also on Antonacci's slides)
% This is the source of all the problems probably. It's way too small for
% PlateZ and since real(X) is what determines the damping the string will
% basically not dampen. For PlateZ2 for some reason its values will break
% everything that comes after this


X = (Y .* 1j .* Z0 .* omega')./pi;

a = zeros(length(f), 300);
b = zeros(length(f), 300);
epsilon = zeros(length(f), 300);

% Detune of the string. Built so that omega*epsilon always goes from -3 to 3 
% as in the graphs
for ii = 1:length(f)
    epsilon(ii,:) = linspace(-0.02*(omega(ii)), 0.02*(omega(ii)), 300);
end

% Acoustics of musical instruments, page 279, eq 6.50 (also on Antonacci's slides)
for ii = 1:length(f)

    a(ii,:) = 1j*imag(X(ii)) + epsilon(ii,:) + sqrt(X(ii).^2 + epsilon(ii,:).^2);
    b(ii,:) = 1j*imag(X(ii)) + epsilon(ii,:) - sqrt(X(ii).^2 + epsilon(ii,:).^2);
    
end


%Plot of eigenfrequency shift wrt string detune. Should be like graphs at
%page 279. For values in PlateZ the graph looks convincing, for PlateZ2 not
%as much: either a or b is way too big or way too small.

figure
for ii = 1:length(f)
    subplot(length(f),1,ii);
    plot(epsilon(ii,:), real(a(ii,:)));
    hold on
    plot(epsilon(ii,:), real(b(ii,:)));
    title(ii + "th freq " + f(ii));
    xlabel("Epsilon [Hz]");
    ylabel("Eigenmode Shift[Hz]");
    legend("a","b");
    grid minor
end

figure
for ii = 1:length(f)
    subplot(length(f),1,ii);
    plot(epsilon(ii,:), imag(a(ii,:)));
    hold on
    plot(epsilon(ii,:), imag(b(ii,:)));
    title(ii + "th freq " + f(ii));
    xlabel("Epsilon [Hz]");
    ylabel("Time Decay");
    legend("a","b");
    grid minor
end

%% Point d

t = linspace(0,1,1000); %time axis

%Displacement computed for an arbitrary epsilon. The higher/lower epsilon
%should give more minima, epsilon values around the middle (150th index)
%should give less minima. For PlateZ if you increase by multiple orders of
%magnitude X the reults are very convincing, for PlateZ2 it seems that the
%damping is way too strong or something


eps_index = 100;

F0 = 1; % F(0), initial condition for force on the strings
mu = sqrt(epsilon(:,eps_index).^2 + (X.^2)); 

X_n = X./(omega');


% paper, eq 19
param_p = X+mu;
param_m = X-mu;

R = abs((param_p.*exp(1j.*param_p.*t) - param_m.*exp(1j.*param_m.*t))./(2.*mu)).^2;

for ii = 1:length(f)
    figure
    plot(t,10*log10(R(ii,:)));
    title(f(ii));
end

%%

tau_a = 1./imag(a);
tau_b = 1./imag(b);
T60 = zeros(length(f),1);
for ii = 1:length(f)

    R_db = 10*log10(R(ii,:)./R(ii,1));
    found = false;
    jj = 1;
    while found == false
        if R_db(jj) <= -60
            T60(ii) = t(jj);
            found = true;
        end
        jj = jj+1;
    end

end
