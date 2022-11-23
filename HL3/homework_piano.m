%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Mirco Pezzoli                                                           %
% 2022                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% Setup
% define the string parameters and the simulation variables defined
% according to the provided values and the numerical implementation.
% We want to implement the finite difference scheme of a piano string tuned
% as C2.

% Temporal sampling parameters
fs = 4*44100;
ts = 1/fs;
signal_length = 8; %[s]
N = signal_length*fs;

% Fundamental note
f_0 = 65.4; %[Hz] C1

% Boundary
zeta_l = 1e20; %[Ohm/(Kg*m^2*s)]
zeta_b = 1000; %[Ohm/(Kg*m^2*s)]

% String parameters
l = 1.92; %[m]
Ms = 35e-3; %[Kg]
Te = 750; %[N]
b1 = 0.25; %[1/s]
b2 = 7.5e-5; %[s]
epsilon = 7.5e-6;

% Spatial sampling parameters
M = 521;

% Aliasing condition


% Number of maximum spatial steps
xm = l/M;

% Integer values
% Spatial sampling
m = 0:1:M;
l_m = l.*m/M;


% FD parameters
courant_number = sqrt(Te/(Ms/l))*(ts/xm);
lambda = courant_number * ts/xm;
c = sqrt(Te/(Ms/l));
mu = epsilon^2/(c^2*xm^2);
nu = (2*b2*ts)/(xm^2);

a1 = (-lambda^2)/(1+b1*ts);
a2 = (lambda^2+4*lambda^2*mu+nu)/(1+b1*ts);
a3 = (2-2*lambda^2-6*lambda^2*mu-2*nu)/(1+b1*ts);
a4 = (-1+b1*ts+2*nu)/(1+b1*ts);
a5 = (-nu)/(1+b1*ts);
af = (ts^2/(Ms/l))/(1+ b1*ts);

% Hammer parameters
Mh = 4.9e-3; %[Kg]
bh = 1e-4; %[s^-1]
d1 = (2)/(1+bh*ts/(2*Mh));
d2 = (-1 + 1+bh*ts/(2*Mh))/(1+bh*ts/(2*Mh));
df = (-ts^2/Mh)/(1+bh*ts/(2*Mh));

% Hammer contact window definition
wh = 0.2;
a = 0.12;
xmh = round(0.12*M);
x0 = a*l;
window_lenght = round(wh*l/xm);
h_hamming_window = hanning(window_lenght); %REMEMBER TO CHECK THIS

h_zeros_a = zeros(round(a*M - window_lenght/2),1);
h_zeros_b = zeros(round(M- a*M - window_lenght/2),1);

h_window = [h_zeros_a' h_hamming_window' h_zeros_b'];

%PDE Coefficients:

% Bridge boundary coefficients
br1 = (2 - 2*lambda^2*mu - 2*lambda^2)/(1 + b1*ts + zeta_b*lambda);
br2 = (4*lambda^2*mu + 2*lambda^2)/(1 + b1*ts + zeta_b*lambda);
br3 = (-2*lambda^2*mu)/(1 + b1*ts + zeta_b*lambda);
br4 = (-1 - b1*ts + zeta_l*lambda)/(1 + b1*ts + zeta_b*lambda);
br5 = (ts^2/(Ms/l))/(1 + b1*ts + zeta_b*lambda);
brf = (ts^2/(Ms/l))/(1+b1*ts+zeta_b*lambda);

% Left hand (hinged string end) boundary coefficients
bl1 = (2 - 2*lambda^2*mu - 2*lambda^2)/(1 + b1*ts + zeta_l*lambda);
bl2 = (4*lambda^2*mu + 2*lambda^2)/(1 + b1*ts + zeta_l*lambda);
bl3 = (-2*lambda^2*mu)/(1 + b1*ts + zeta_l*lambda);
bl4 = (-1 - b1*ts + zeta_l*lambda)/(1 + b1*ts + zeta_l*lambda);
bl5 = (ts^2/(Ms/l))/(1 + b1*ts + zeta_l*lambda);
blf = (ts^2/(Ms/l))/(1+b1*ts+zeta_l*lambda);

% Hammer felt parameters
ph = 2.3; %Felt stiffness exponent
Kh = 4e8; %Felt stiffness 
%epsilon = |hammer displacement - string displacement at hammer
%point|(EQ 3.27)

%% Computation of the FD scheme
% Initialization
y = zeros(M,N);
eta = zeros(1,N);
Fh = zeros(1,N);
F = zeros(M,N);

for nn = 1:N
    Fh(nn) = Kh*abs(eta - y(xmh,nn))^ph;
    eta(0,nn+1) = d1*eta(0,nn)+d2*eta(0,nn-1)+df *Fh(nn);

    F(:,nn) = Fh.*h_window;

    y(0,nn+1) = bl1.*y(0,nn) + ...
                bl2.*y(1,nn) + ...
                bl3.*y(2,nn) + ...
                bl4.*y(0,nn-1) + ...
                blf*Fh(nn);

    y(1,nn+1) = a1*(y(3,nn) - y(1,nn) + 2*y(0,nn)) + ...
                a2*(y(2,nn) + y(0,nn)) + ...
                a3*y(1,nn) + ...
                a4*y(1,nn-1) + ...
                a5*(y(2,nn-1) + y(0,nn-1)) + ...
                af * F(1,nn);
    for mm = 2:M-2
        
        y(mm,nn+1) = a1*(y(mm+2,nn) + y(mm-2,nn)) + ...
                     a2*(y(mm+1,nn) + y(mm-1,nn)) + ...
                     a3*y(mm,nn) + ...
                     a4*y(mm,nn-1) + ...
                     a5*(y(mm+1,nn-1)+y(mm-1,nn-1)) + ...
                     af*F(mm,nn);
    end

    y(M-1,nn) = a1*(2*y(M,nn) - y(M-1, nn) + y(M-3, nn)) + ...
                a2*(y(M,nn) + y(M-2,nn)) + ...
                a3*(y(M-1, nn)) + ...
                a4*y(M-1,nn-1) + ...
                a5*(y(M,nn-1) + y(M-2,nn-1)) + ...
                af*F(M-1,nn);

    y(M,nn) = br1*y(M,nn) + br2*y(M-1, nn) + br2*y(M-2,nn) + br4*y(M,nn-1) + brf*F(M,nn);x

end


% Computation loop
%% Plot the displacement in time

%% Plot the synthesized signal play it and save it on the disk

% Play the sound



% Save on disk










