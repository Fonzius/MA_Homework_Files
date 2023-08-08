clear all;
clc;
close all;
%% Data
r = 0.15; %[m]
T = 10; %[N/m]
sigma = 0.07; %[kg/m^2]

%% part 1, point a

c = sqrt(T/sigma); %[m/s]

%% part 1, point b

bessel_zeros = [2.4048, 3.8317, 5.1356 , 5.5201, 6.3802, 7.0156, ...
                7.5883, 8.4172, 8.6537, 8.7715 , 9.7610, 9.9361, ...
                10.1735, 11.0647, 11.0864, 11.6198, 11.7915, 12.2251];
% REMEMBER TO ASSUME THAT THE MEMBRANE IS CLAMPED

w = (bessel_zeros./(r)).*c;

m = [0, 1, 2, 0, 3, 1, 4, 2, 0, 1, 5, 3, 6, 1, 4, 7, 2, 8]; %Mode numebrs (from the table)
n = [1, 1, 1, 2, 1, 2, 1, 2, 3, 6, 1, 2, 1, 3, 2, 1, 3, 1];
phi = linspace(0, 2*pi, 1000);
radius = linspace(0,r, 1000);

k_a = w./c;

for ii = 1:6
    z_part1 = exp(-1i .* m(ii) .* phi);
    z_part2 = besselj(m(ii), k_a(ii).*radius);
    
    z = real(z_part1 .* z_part2');
    
    %figure()
    subplot(2,3,ii);
    polarplot3d(z, 'RadialRange', [0 r], 'TickSpacing', 0);
    title("Modeshape_" + m(ii)+"_"+n(ii));
    colorbar;
    %shading interp % this changes the color shading (i.e. gets rid of the grids lines on the surface)
   view(2) % same as view(0, 90)
end

sgtitle("First 6 modes of vibration of the circular membrane")

%%
r_in = 0.075; %[m]
r_out = 0.075; %[m]
phi_in = 15; %[deg]
phi_out = 195; %[deg]

%Translate angles to radiants
%phi_in = phi_in * 2*pi / 360;
%phi_out = phi_out * 2*pi / 360;

Q_modes = 25;

time = linspace(0,2,10000);
input_force = 0.1.*exp(-((time-0.03).^2)./0.001^2);

eigen_angular_freqs = (bessel_zeros./r)*c;

omega = linspace(min(eigen_angular_freqs)/10, max(eigen_angular_freqs)*2, 10000);

modal_matrix=zeros(18,10000);
z_in = zeros(18,1);
z_out = zeros(18,1);

mass = r^2*pi*sigma;

for ii = 1:18
    if m(ii) ~= 0
        phi_in = pi/m(ii);
        phi_out = phi_in + pi;
    else
        phi_in = 0;
        phi_out = 180;
    end
    modal_matrix(ii,:) = abs(1./(-omega.^2 + 1i.*eigen_angular_freqs(ii).*omega./Q_modes + eigen_angular_freqs(ii)^2));
    z_in(ii) = exp(-1i * m(ii) * phi_in) * besselj(m(ii), k_a(ii)*r_in);
    z_out(ii) = exp(-1i * m(ii) * phi_out) * besselj(m(ii), k_a(ii)*r_out);
end

output_modes = zeros(10000,1);

for ii = 1:10000
    modal_temp = diag(modal_matrix(:, ii));
    output_modes(ii) = z_in' * modal_temp * z_out;
end


figure()
plot(omega,abs(output_modes));

force_fft = fft(input_force,10000);
output_fft = force_fft .* output_modes;

output_displacement = ifft(output_fft);
figure()
plot(time(1:300), output_displacement(1:300));



%% Part 2, point d
h = 0.001; %[m]
E = 69*10^9; %[pa]
rho = 2700; %[kg/m^3]
mu = 0.334;

%quasi-longitudinal wave speed
cl = sqrt(E/(rho*(1-mu^2)));
%longitudinal waves speed
cl_ = sqrt((E*(1-mu))/(rho*(1+mu)*(1-2*mu)));

%% Part 2, point e
freq_axis = 20:1:20000;
v = sqrt(1.8*freq_axis*h*cl);

figure();
semilogx(freq_axis,v, 'b','LineWidth',2);
title("Propagation speed of bending waves");
ylabel("$\rm{Speed\quad}$[m/s]",'interpreter', 'latex');
xlabel("$\rm{Frequency\quad}$[Hz]", 'interpreter', 'latex');

%% Part 2, point f
%plate bending modes
f_01 = (0.4694*cl*h)/(r^2);
f_11 = 2.08*f_01;
f_21 = 3.41*f_01;
f_02 = 3.89*f_01;
f_31 = 5.00*f_01;

f_0s_plate = [f_01, f_11, f_21, f_02, f_31];

%% Part 3, point g

rho_string = 5000; %[kg/m^3]
r_string = 0.001; %[m]
l_string = 0.4; %[m]

Q = 50;

T_string = rho_string*r_string^2*pi*((f_01*2*pi*l_string)^2)/(1^2*pi^2);


%% Part 3, point h

E_string = 200*10^9;

K_string = r_string/2;
S_string = r_string^2*pi;

B = pi*E_string*S_string*K_string^2/(T_string*l_string^2);

f_0_string = f_01;

% assuming supported end
f_0_stiff_string = zeros(5,1);
for ii = 1:5
    %f_0_stiff_string(ii) = ii*f_0_string*sqrt(1+B*ii^2);
    f_0_stiff_string(ii) = ii*f_0_string*sqrt(1+B*ii^2)*(1 + 2*sqrt(B)/pi + 4*B/pi^2);
end

%% Part 3, point i

m_string = S_string*l_string* rho_string;
m_plate = rho*h*r^2*pi;

modes =  [1:1:10];

coupling1 = m_string./(modes.^2*m_plate);
coupling2 = (pi^2) / (4*Q^2);

coupling = coupling1 - coupling2; % if > 0 the modes are strogly coupled, if < 0 weakly coupled

f_0s_string = zeros(5,1); % first 5 natural frequencies of the string (stiffness is not kept into consideration)
for ii = 1:5
    f_0s_string(ii) = ii * f_0_string; 
end

deltaX = zeros(5:1);
for ii = 1:5
    deltaX(ii) = (f_0s_string(ii)-f_0s_plate(ii))/f_0s_plate(ii); 
end

% computed with processing
deltaF2 = [0.009999998 , -0.046000004];
deltaF4 = [0.0388 , -0.010399997];

f_2_coupled_string = deltaF2.*f_0s_plate(2)+f_0s_plate(2);
f_4_coupled_string = deltaF4.*f_0s_plate(4)+f_0s_plate(4);

coupling_coefficient = coupling1.*10^3;

% computed, once again, with processing
f_1_coupled = [f_0s_string(1) + 0.0105999*f_0s_string(1)/2, f_0s_string(1) - 0.0105999*f_0s_string(1)/2];


