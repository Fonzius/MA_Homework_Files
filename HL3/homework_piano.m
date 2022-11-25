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
rho = Ms/l;
Te = 4*l^2*rho*f_0^2;%750; %[N]
b1 = 0.5; %[1/s]
b2 = 6.25e-9; %[s]
epsilon = 7.5e-6;
c = sqrt(Te/rho);


% Spatial sampling parameters
M = 520;

% Aliasing condition
F_ny = fs/2;
gamma = F_ny/f_0;

M_max = l / sqrt(0.5*(c^2*ts^2+4*b2*ts+sqrt((c^2*ts^2+4*b2*ts)^2+16*epsilon^2*ts^2))); 
M_min = ((-1+(1+16*epsilon*gamma^2)^(1/2))/(8*epsilon))^(1/2);

% M_min < M=520 < M_max, no aliasing

% Number of maximum spatial steps
xm = l/M;

% Integer values
% Spatial sampling
m = 0:1:M;
l_m = l.*m/M;


% FD parameters
%courant_number = c*(ts/xm);
lambda = c * ts/xm;
mu = epsilon^2/(c^2*xm^2);
nu = (2*b2*ts)/(xm^2);

a1 = (-lambda^2*mu)/(1+b1*ts);
a2 = (lambda^2 + 4*lambda^2*mu + nu)/(1 + b1*ts);
a3 = (2 - 2*lambda^2 - 6*lambda^2*mu - 2*nu)/(1+b1*ts);
a4 = (-1 + b1*ts + 2*nu)/(1+b1*ts);
a5 = (-nu)/(1+b1*ts);
af = (ts^2/rho)/(1+ b1*ts);

% Hammer parameters
Mh = 4.9e-3; %[Kg]
bh = 1e-4; %[s^-1]
d1 = (2)/(1+bh*ts/(2*Mh));
d2 = (-1 + bh*ts/(2*Mh))/(1+bh*ts/(2*Mh));
df = (-ts^2/Mh)/(1+bh*ts/(2*Mh));
v0h = 2.5; %[m/s]

% Hammer contact window definition
wh = 0.2;
a = 0.12;
xmh = round(0.12*M);
x0 = a*l;
window_lenght = round(wh/xm);
h_window = zeros(1,M);
h_window(1, xmh-window_lenght/2 : xmh+window_lenght/2-1) = hanning(window_lenght);

%PDE Coefficients:

% Bridge boundary coefficients
br1 = (2 - 2*lambda^2*mu - 2*lambda^2)/(1 + b1*ts + zeta_b*lambda);
br2 = (4*lambda^2*mu + 2*lambda^2)/(1 + b1*ts + zeta_b*lambda);
br3 = (-2*lambda^2*mu)/(1 + b1*ts + zeta_b*lambda);
br4 = (-1 + b1*ts + zeta_b*lambda)/(1 + b1*ts + zeta_b*lambda);
brf = (ts^2/rho)/(1 + b1*ts + zeta_b*lambda);

% Left hand (hinged string end) boundary coefficients
bl1 = (2 - 2*lambda^2*mu - 2*lambda^2)/(1 + b1*ts + zeta_l*lambda);
bl2 = (4*lambda^2*mu + 2*lambda^2)/(1 + b1*ts + zeta_l*lambda);
bl3 = (-2*lambda^2*mu)/(1 + b1*ts + zeta_l*lambda);
bl4 = (-1 - b1*ts + zeta_l*lambda)/(1 + b1*ts + zeta_l*lambda);
blf = (ts^2/rho)/(1 + b1*ts + zeta_l*lambda);

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

avg_disp = zeros(1,N);

%Step 2 ( t=1 )
eta(1,2) = v0h*ts;
Fh(1,2) = Kh.*abs(eta(1,1) - y(xmh,1)).^ph;

%Step 3 (t=2)
for mm = 2:M-1
    y(mm,3) = y(mm-1,2)+y(mm+1,2)-y(mm,1)+(ts^2.*M.*Fh(1,2).*h_window(1,mm))/Ms;% see Chaigne pg.5
end
eta(1,4) = 2*eta(1,2) - eta(1,1) - (ts^2.*Fh(1,2))/Mh;
eta(1,3) = eta(1,4);


% Computation loop

for nn = 4:N
    
    if eta(1,nn) < y(xmh,nn)
        Fh(1,nn) = 0;
    else
        Fh(1,nn) = Kh.*abs(eta(1,nn) - y(xmh,nn)).^ph;
    end
    
    F(:,nn) = Fh(1,nn).*h_window(1,:);
    

    y(1,nn+1) = bl1.*y(1,nn) + ...
                bl2.*y(2,nn) + ...
                bl3.*y(3,nn) + ...
                bl4.*y(1,nn-1) + ...
                blf.*F(1,nn);

    y(2,nn+1) = a1.*(y(4,nn) - y(2,nn) + 2.*y(1,nn)) + ...
                a2.*(y(3,nn) + y(1,nn)) + ...
                a3.*y(2,nn) + ...
                a4.*y(2,nn-1) + ...
                a5.*(y(3,nn-1) + y(1,nn-1)) + ...
                af.*F(1,nn);

    for mm = 3:M-2
        
        y(mm,nn+1) = a1.*(y(mm+2,nn) + y(mm-2,nn)) + ...
                     a2.*(y(mm+1,nn) + y(mm-1,nn)) + ...
                     a3.*y(mm,nn) + ...
                     a4.*y(mm,nn-1) + ...
                     a5.*(y(mm+1,nn-1)+y(mm-1,nn-1)) + ...
                     af.*F(mm,nn);
    end

    y(M-1,nn) = a1.*(2*y(M,nn) - y(M-1, nn) + y(M-3, nn)) + ...
                a2.*(y(M,nn) + y(M-2,nn)) + ...
                a3.*(y(M-1, nn)) + ...
                a4.*y(M-1,nn-1) + ...
                a5.*(y(M,nn-1) + y(M-2,nn-1)) + ...
                af.*F(M-1,nn);

    y(M,nn) = br1.*y(M,nn) + ...
              br2.*y(M-1, nn) + ...
              br2.*y(M-2,nn) + ...
              br4.*y(M,nn-1) + ...
              brf.*F(M,nn);

    eta(1,nn+1) = d1.*eta(1,nn) + d2.*eta(1,nn-1) + df.*Fh(1,nn);
    
    for mavg = M-12:M
        avg_disp(1,nn) = avg_disp(1,nn) + y(mavg, nn);
    end
    
    avg_disp(1,nn) = avg_disp(1,nn)/12;
end

%% Plot the displacement in time

figure
for ii = 1:100:size(y,2)
    plot(y(:,ii),'LineWidth',2)
    ylim([-0.0000001 0.0000001])
    xlim([0,M]);
    xlabel('x');
    ylabel("y");
    title(ii);
    pause(0.000001);
end

%% Plot the synthesized signal, play it and save it on the disk

% Play the sound

t = 0:ts:signal_length-ts;
figure
plot(t, avg_disp, LineWidth=2);

soundsc(avg_disp, fs);
% Save on disk
sound = decimate(avg_disp, 4, 'fir');
sound = (sound/max(abs(sound)));
audiowrite('10669941_Bernasconi_10876787_Luan_Piano.wav', sound, fs/4);









