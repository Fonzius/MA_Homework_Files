clc;
clear all;
close all;

%% Point b
close all;

L = 1;
H = 1.4;

ccc = [0 0.4470 0.7410]; %plot color

f = [349.23, 440, 523.25, 659.25, 783.99];

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
subplot(5,1,1)
plot(Z_f1(1:10,1),Z_f1(1:10,2), '-o', Color='r');
hold on
plot(L-Z_f1(1:10,1),Z_f1(1:10,2), '-o', Color='r');
hold on
plot(Z_f1(1:10,1),H-Z_f1(1:10,2), '-o', Color='r');
hold on
plot(L-Z_f1(1:10,1),H-Z_f1(1:10,2), '-o', Color='r');
title(f(1) + " [Hz]");
xlabel("x[m]");
ylabel("y[m]");

subplot(5,1,2);
plot(Z_f2(1:10,1),Z_f2(1:10,2), '-o', Color='g');
hold on
plot(L-Z_f2(1:10,1),Z_f2(1:10,2), '-o', Color='g');
hold on
plot(Z_f2(1:10,1),H-Z_f2(1:10,2), '-o', Color='g');
hold on
plot(L-Z_f2(1:10,1),H-Z_f2(1:10,2), '-o', Color='g');
title(f(2) + " [Hz]");
xlabel("x[m]");
ylabel("y[m]");

subplot(5,1,3);
plot(Z_f3(1:10,1),Z_f3(1:10,2), '-o', Color='b');
hold on
plot(L-Z_f3(1:10,1),Z_f3(1:10,2), '-o', Color='b');
hold on
plot(Z_f3(1:10,1),H-Z_f3(1:10,2), '-o', Color='b');
hold on
plot(L-Z_f3(1:10,1),H-Z_f3(1:10,2), '-o', Color='b');
title(f(3) + " [Hz]");
xlabel("x[m]");
ylabel("y[m]");

subplot(5,1,4);
plot(Z_f4(1:10,1),Z_f4(1:10,2), '-o', Color='k');
hold on
plot(L-Z_f4(1:10,1),Z_f4(1:10,2), '-o', Color='k');
hold on
plot(Z_f4(1:10,1),H-Z_f4(1:10,2), '-o', Color='k');
hold on
plot(L-Z_f4(1:10,1),H-Z_f4(1:10,2), '-o', Color='k');
title(f(4) + " [Hz]");
xlabel("x[m]");
ylabel("y[m]");

subplot(5,1,5);
plot(Z_f5(10:20,1),Z_f5(1:10,2), '-o', Color='magenta');
hold on
plot(L-Z_f5(1:10,1),Z_f5(1:10,2), '-o', Color='magenta');
hold on
plot(Z_f5(1:10,1),H-Z_f5(1:10,2), '-o', Color='magenta');
hold on
plot(L-Z_f5(1:10,1),H-Z_f5(1:10,2), '-o', Color='magenta');
title(f(5) + " [Hz]");
xlabel("x[m]");
ylabel("y[m]");


Z_final = [L-Z_f1(7,1), H-Z_f1(7,2), Z_f1(7,3:5); ...
           L-Z_f2(4,1), H-Z_f2(4,2), Z_f2(4,3:5); ...
           L-Z_f3(10,1), Z_f3(10,2), Z_f3(10,3:5); ...
           L-Z_f4(4,1), Z_f4(4,2), Z_f4(4,3:5); ...
           Z_f5(7,1), Z_f5(7,2), Z_f5(7,3:5); ...
           ];
figure
plot(Z_final(:,1), Z_final(:,2), '-o', LineWidth=2);
axis([0,L, 0,H]);
xlabel("x[m]");
ylabel("y[m]");
%title("Bridge Design");

%% Point c

rho = 10.8e-3; %[kg/m] piano wire linear density, from the HW text

f = [349.23, 440, 523.25, 659.25, 783.99];
omega = f.*2.*pi;

T = [800;800;800;800;800];
Z0 = sqrt(T.*rho);
Y = 1./(Z_final(:,5));

% Acoustics of musical instruments, page 275, eq 6.34 (also on Antonacci's slides)

X = (Y .* 1j .* Z0 .* omega')./pi;

a = zeros(length(f), 300);
b = zeros(length(f), 300);
epsilon = zeros(length(f), 300);


for ii = 1:length(f)
    epsilon(ii,:) = linspace(-0.02*(omega(ii)), 0.02*(omega(ii)), 300);
end

% Acoustics of musical instruments, page 279, eq 6.50 (also on Antonacci's slides)
for ii = 1:length(f)

    a(ii,:) = 1j*imag(X(ii)) + epsilon(ii,:) + sqrt(X(ii).^2 + epsilon(ii,:).^2);
    b(ii,:) = 1j*imag(X(ii)) + epsilon(ii,:) - sqrt(X(ii).^2 + epsilon(ii,:).^2);
    
end

figure
for ii = 1:length(f)
    subplot(length(f),1,ii);
    plot(epsilon(ii,:), real(a(ii,:)), Color=ccc);
    hold on
    plot(epsilon(ii,:), real(b(ii,:)), Color=ccc);
    title(f(ii) + " [Hz]");
    xlabel("$\epsilon$ [Hz]", Interpreter="latex");
    ylabel("$\Re{[a]}$[Hz]", Interpreter="latex");
    xlim([epsilon(ii,1) epsilon(ii,300)]);
    %legend("a","b");
    grid minor
end

figure
for ii = 1:length(f)
    subplot(length(f),1,ii);
    plot(epsilon(ii,:), imag(a(ii,:)), Color=ccc);
    hold on
    plot(epsilon(ii,:), imag(b(ii,:)), Color=ccc);
    title(f(ii) + " [Hz]");
    xlabel("$\epsilon$ [Hz]", Interpreter="latex");
    ylabel("$\Im{[a]}$[$s^{-1}$]", Interpreter="latex");
    xlim([epsilon(ii,1) epsilon(ii,300)]);
    %legend("a","b");
    grid minor
end

%% Point d
close all;

t_int = 60000;

t = linspace(0,80,t_int); %time axis

%eps_index = [150, 200, 250];

R = zeros(5,t_int,3);

for ii = 1:300
    mu = sqrt(epsilon(:,ii).^2 + (X.^2)); 
    
    X_n = X./(omega');
    
    % paper, eq 19
    param_p = X+mu;
    param_m = X-mu;
    
    R(:,:,ii) = abs((param_p.*exp(1j.*param_p.*t) - param_m.*exp(1j.*param_m.*t))./(2.*mu)).^2;
end

figure
position = 1;

for ii = 1:5
    for jj = 150:50:300
        subplot(5,4,position)
        plot(t,10*log10(R(ii,:,jj)));
        title(f(ii) + " [Hz], $\epsilon$ = " + epsilon(ii, jj) + " [Hz]", Interpreter="latex");
        xlim([0 1]);

        position = position+1;
    end
end

T60 = zeros(size(epsilon));

for ii = 1:5
    for jj = 1:300

        R_db = 10*log10(R(ii,:,jj)./R(ii,1,jj));
        found = false;
        kk = 1;
        while found == false && kk <= t_int
            if R_db(kk) <= -60
                T60(ii,jj) = t(kk);
                found = true;
            end
            kk = kk+1;
        end
    end
end

figure
for ii = 1:5
    subplot(5,1,ii)
    plot(epsilon(ii,:), T60(ii,:));
    xlim([epsilon(ii,1) epsilon(ii,300)]);
    title(f(ii) + " [Hz]");
    grid minor;
end


%% TIME CONSTANTS
close all;

tau_a = 1./imag(a);
tau_b = 1./imag(b);

figure
subplot(2,1,1)
semilogy(epsilon,tau_a);
xlabel("$\epsilon$ [Hz]", Interpreter="latex");
ylabel("$\tau_+$[s]]", Interpreter="latex");
subplot(2,1,2)
semilogy(epsilon, tau_b);
xlabel("$\epsilon$ [Hz]", Interpreter="latex");
ylabel("$\tau_-$[s]", Interpreter="latex");

%% FULL FREQUENCY PLOTS
close all;

frm = importdata("FRM.txt");

for ii=1:3 % Removing empty rows
    frm(1,:) = [];
end

%frm_sorted = zeros(35,251);
frm_f = str2double(frm(1:251,3));       %full frequency axis

points = frm(1:251:8785, 1:2);                        %list of points
frm_sorted = reshape(str2double(frm(:,4)),[251,35])'; %each row has the frf at one point

figure
for ii = 1:length(points(:,1))
    subplot(5,7,ii);
    plot(frm_f,frm_sorted(ii,:));
    xlim([frm_f(1) frm_f(251)]);
    ylim([30,80]);
    xlabel("x: " + points(ii,1) + " y: " + points(ii,2));
    
end





