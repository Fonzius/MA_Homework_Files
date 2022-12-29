clc;
clear all;
close all;
%% Question 1
alpha = 0.75 * pi/180; %[rad]
L = 0.45; %[m]
f_e4 = 329.63; %[hz]
c = 343; %[m/s]
rho = 1.225;

L_c = 0.04; %[m] Mouth end correction for low freqs

x1 = ((pi*c)/(2*pi*f_e4*(1+0.85*tan(alpha))))*(L+L_c)*(1+0.85*tan(alpha)); %Derived from wn = (n*pi*c)/(L+x1)
a1 = x1*tan(alpha); %Mouth radius
a2 = (x1+L)*tan(alpha); %End radius

%% Question 1, Impedance
r1 = linspace(L*tan(alpha),0.01, 1000);
r0 = r1 -L*tan(alpha);
Z_full = zeros(1,1000);
k_e4 = 2*pi*f_e4/c;

for ii = 1:length(r1)
    Z_L0=0.25.*(f_e4*2*pi).^2.*rho./pi./c+1j.*0.61.*rho.*f_e4*2*pi/pi./r1(ii);
    Z_full(ii) = ZIN1(r0(ii), r1(ii),L+0.61*L+L_c, Z_L0,  k_e4, rho, c);
end

[Z_full_max, max_full_index] = max(abs(Z_full));
r0_max = r0(max_full_index);
r1_max = r1(max_full_index);


%% Question 2

f_fs4 = 369.99; %[hz]

k_fs4 = 2*pi*f_fs4/c;
S1 = pi*a1^2;
S2 = pi*a2^2;

x2 = x1+L;

D = linspace(a2,L-a2,200);

Z_total = zeros(1,length(D));

Zl = 1j*k_fs4*rho*c*0.85*a2/S2;
Z_end = zeros(1,length(D));
Z_hole = zeros(1,length(D));
Z_ahah = zeros(1,length(D));

for ii = 1:length(D)
    r1 = (L-D(ii))*tan(alpha);

    Z_end(ii) = ZIN1(r1, a2, D(ii), Zl, k_fs4, rho, c);
    Z_hole(ii) = 1/(-1j* a2*cot(k_fs4*0.85*L/(rho*c))); %Acoustic length of the hole simplified to be 0.85L as page 466 Rossing
    Z_ahah(ii) = 1./((1./Z_end(ii)) + (1./Z_hole(ii)));

    Z_total(ii) = ZIN1(a1, r1, L-D(ii), Z_ahah(ii), k_fs4, rho, c);
end

[Z_max, max_index] = max(abs(Z_total));
D_max = D(max_index);

%% Question 3
f_gs4 = 415.3; %[Hz]
k_gs4 = 2*pi*f_gs4/c;

Z_in_g = Z_ahah(max_index);
D_g = linspace(D_max+2*a2,L-a2,200);
r_2 = (L-D_max)*tan(alpha);

Z_end_g = ZIN1(r1, a2, D(ii), Zl, k_gs4, rho, c);
Z_hole_g = 1/(-1j* a2*cot(k_gs4*0.85*L/(rho*c)));
Z_ending_g = 1./((1./Z_end_g) + (1./Z_hole_g));

Z_total_g = zeros(1, length(D_g));

for ii = 1:length(D)
    r1 = (L-D(ii))*tan(alpha);
    
    Z_middle_g = ZIN1(r1, r_2, D_g(ii)-D_max, Z_ending_g, k_gs4, rho,c);
    Z_total_middle_g = 1./((1./Z_middle_g) + (1./Z_hole_g));
    Z_total_g(ii) = ZIN1(a1, r1, L-D_g(ii), Z_total_middle_g, k_gs4, rho,c);
end

[Z_max_g, max_index_g] = max(abs(Z_total_g));
D_max_g = D_g(max_index_g);

%check cut frequency, worst case pipe with fixed radius a2

f_c = 0.11*c*sqrt(1/(D_max_g-D_max));

%% Question 4
delta_p = 55; %[Pa]
f_centroid = 1700; %[Hz]
nu = 1.5e-5; %Kinematic viscosity (jetswitch instrument pt1, slide 21)

Uj = sqrt(2*delta_p/rho);

h = 0.3*Uj/f_centroid;

Re = Uj*h/nu;

%for values of the Reynolds number smaller than 2000, the jet remains laminar for a short distance;
%Jetswitch instruments pt2, page 2

%% Question 5
L_channel = 0.025; %[m]
channel_thickness = sqrt(nu*L_channel/Uj);

