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
% for ii=1:3 %Removing not needed data
%     Z(1,:) = [];
% end

Z = importdata("PlateZ2.txt");

for ii=1:4 %Removing not needed data
    Z(1,:) = [];
end

% 1 = x
% 2 = y
% 3 = f
% 4 = abs(Z)
% 5 = Z

Z_sorted = sortrows(Z,3);

%Sort for each frequency by db(abs(Z))
Z_f1 = str2double(sortrows(Z_sorted(1:48,:),4));
Z_f2 = str2double(sortrows(Z_sorted(49:96,:),4));
Z_f3 = str2double(sortrows(Z_sorted(97:144,:),4));
Z_f4 = str2double(sortrows(Z_sorted(145:192,:),4));
Z_f5 = str2double(sortrows(Z_sorted(193:240,:),4));

%ii_f = ones(5,1);

%found_best = 0;
% 
% while found_best == 0
%     found_best = 1;
%     Z_best = [Z_f1(ii_f(1),:); ...
%               Z_f2(ii_f(2),:); ...
%               Z_f3(ii_f(3),:); ...
%               Z_f4(ii_f(4),:); ...
%               Z_f5(ii_f(5),:)];
%     
%     Z_corrected = str2double(Z_best);
%     
%     for ii = 2:5
%         if Z_corrected(ii,1) < Z_corrected(ii-1,1)
%             Z_corrected(ii,1) = Z_corrected(ii,1) + L/2;
%         end
%     end
%     
%     for ii = 2:5
%         if Z_corrected(ii,1) < Z_corrected(ii-1,1)
%             ii_f(ii-1,1) = ii_f(ii-1,1) + 1;
%             found_best = 0;
%         end
%     end
% 
% end

figure
plot(Z_f1(1:6,1),H-Z_f1(1:6,2), '-o', Color='r');
hold on
plot(Z_f2(1:6,1),H-Z_f2(1:6,2), '-o', Color='g');
hold on
plot(Z_f3(1:6,1),H-Z_f3(1:6,2), '-o', Color='b');
hold on
plot(Z_f4(1:6,1),H-Z_f4(1:6,2), '-o', Color='k');
hold on
plot(Z_f5(1:7,1),H-Z_f5(1:7,2), '-o', Color='magenta');
hold on

plot(L-Z_f1(1:6,1),H-Z_f1(1:6,2), '-o', Color='r');
hold on
plot(L-Z_f2(1:6,1),H-Z_f2(1:6,2), '-o', Color='g');
hold on
plot(L-Z_f3(1:6,1),H-Z_f3(1:6,2), '-o', Color='b');
hold on
plot(L-Z_f4(1:6,1),H-Z_f4(1:6,2), '-o', Color='k');
hold on
plot(L-Z_f5(1:7,1),H-Z_f5(1:7,2), '-o', Color='magenta');


% figure
% plot(Z_f1(1,1),Z_f1(1,2)+H/2, '-o', Color='r');
% hold on
% plot(Z_f2(5,1),Z_f2(5,2)+H/2, '-o', Color='g');
% hold on
% plot(Z_f3(6,1),Z_f3(6,2)+H/2, '-o', Color='b');
% hold on
% plot(Z_f4(9,1)+L/2,Z_f4(9,2), '-o', Color='k');
% hold on
% plot(Z_f5(5,1)+L/2,Z_f5(5,2), '-o', Color='magenta');
% axis([0,L, 0,H]);

Z_final = [Z_f1(2,1), H-Z_f1(2,2), Z_f1(2,3:5); ...
           Z_f2(1,1), H-Z_f2(1,2), Z_f2(1,3:5); ...
           Z_f3(6,1), H-Z_f3(6,2), Z_f3(6,3:5); ...
           L-Z_f4(2,1), Z_f4(2,2), Z_f4(2,3:5); ...
           L-Z_f5(6,1), Z_f5(6,1), Z_f5(6,3:5); ...
           ];
figure
plot(Z_final(:,1), Z_final(:,2));
axis([0,L, 0,H]);

%% Point c
rho = 0.0108; %[kg/m]

f = [349.23, 440, 523.25, 659.25, 783.99];
omega = f.*2.*pi;

c = 2.*Z_final(:,2).*f';
T = rho.*c.^2;

Z0 = T./c;
Y = 1./(Z_final(:,5));

X = Y .* (1j .* Z0)./pi;

a = zeros(length(f), 300);
b = zeros(length(f), 300);
epsilon = zeros(length(f), 300);

for ii = 1:length(f)
    epsilon(ii,:) = linspace(-3/(omega(ii)), 3/(omega(ii)), 300);
end

for ii = 1:length(f)

    a(ii,:) = X(ii) + epsilon(ii,:) + sqrt(X(ii).^2 + epsilon(ii,:).^2);
    b(ii,:) = X(ii) + epsilon(ii,:) - sqrt(X(ii).^2 + epsilon(ii,:).^2);
    
end

figure
for ii = 1:length(f)
    subplot(length(f),1,ii);
    plot(omega(ii)*epsilon(ii,:), omega(ii)*a(ii,:));
    hold on
    plot(omega(ii)*epsilon(ii,:), omega(ii)*b(ii,:));
    title(ii + "th freq " + f(ii));
end

%% Point d

t = linspace(0,10,1000);

x = zeros(length(f), 1000);
for ii = 1:length(f)
    x(ii,:) = linspace(0,Z_final(ii,2),1000);
end

X = X.*1000000;
F0 = 1;
mu = sqrt(epsilon(:,10).^2 + (X.^2));

Vb = ((2.*pi.*F0.*X)./(mu.*Z0)) .* exp(1j.*(epsilon(:,10) + X + 1) .* (omega'.*t)) .* (mu.*cos(mu.*omega'.*t) + 1j.*X.*sin(mu.*omega'.*t));

figure
plot(t, db(abs(Vb(5,:)./Vb(5,1))));

