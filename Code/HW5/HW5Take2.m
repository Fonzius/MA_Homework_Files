clc;
clear all;
close all;

alpha = 0.75 * pi/180; %[rad]
L = 0.45; %[m]
f_e4 = 329.63; %[hz]
c = 343; %[m/s]
rho = 1.225;

L_c = 0.04; %[m] Mouth end correction for low freqs

%% Question 1
close all;

r2 = linspace(L*tan(alpha),0.1, 10000);
r1 = r2 - L*tan(alpha);

k_e4 = 2*pi*f_e4/c;
omega_e4 = 2*pi*f_e4;

Z_mouth_e4 = 1j* (L_c*rho./(pi.*r2.^2)) .* 2.*pi.*f_e4;

%Z_L_e4 = (0.25.*(f_e4*2*pi).^2.*rho./pi./c+1j.*0.61.*rho.*f_e4*2*pi/pi./r1).*((1+cos(alpha))/2);

Z_e4 = ZIN1(r2, r1, L + 0.85.*r1, 0,  k_e4, rho, c) + Z_mouth_e4;

max_index_e4 = find(islocalmin(db(Z_e4)),1,'last');
a2 = r2(max_index_e4); %radius at the mouth of the player
a1 = r1(max_index_e4); %radius at the end of the recorder
Z_max = Z_e4(max_index_e4);
x1 = a2/tan(alpha);

L_corr = L+a1*0.85;

figure
plot(r1,db(Z_e4));
xline(r1(max_index_e4));
xlabel("Foot Radius[m]");
xlim([r1(1), r1(10000)]);
ylabel("Z dB");

%% Question 2
close all;

f_fs4 = 369.99; %[hz]

k_fs4 = 2*pi*f_fs4/c;

Z_mouth_fs4 = 1j* (L_c*rho./(pi.*a2.^2)) .* 2.*pi.*f_fs4;

%Z_L_fs4= (0.25.*(f_fs4*2*pi).^2.*rho./pi./c+1j.*0.61.*rho.*f_fs4*2*pi/pi./a1).*((1+cos(alpha))/2);

D_fs4 = linspace(a1,L-a1,10000); %goes from mouth of the player to foot of the recorder

r_fs4 = a2 - D_fs4.*tan(alpha); %radius of the pipe where the hole is

Z_hole_fs4 = 1./(-1j* ((pi*a1^2)/(rho*c))*cot(k_fs4.*(0.85.*a1+r_fs4))); %Acoustic length of the hole simplified to be 0.85L as page 466 Rossing

Z_in_fs4 = ZIN1(r_fs4, a1, L_corr-D_fs4, 0,  k_fs4, rho, c);
Z_end_fs4 = 1./((1./(Z_in_fs4)) + (1./(Z_hole_fs4)));
Z_fs4 = ZIN1(a2, r_fs4, D_fs4, Z_end_fs4,  k_fs4, rho, c) + Z_mouth_fs4;

max_index_fs4 = find(islocalmin(db(Z_fs4)), 1, 'last');

x_fs4 = D_fs4(max_index_fs4);
r_max_fs4 = r_fs4(max_index_fs4);

figure
plot(D_fs4,db(Z_fs4));
xline(D_fs4(max_index_fs4));
xlabel("x[m]");
xlim([D_fs4(1), D_fs4(10000)]);
ylabel("Z dB");


%% Question 3

f_gs4 = 415.3; %[Hz]
k_gs4 = 2*pi*f_gs4/c;

D_gs4 = linspace(a1,x_fs4-2*a1,10000);
r_gs4 = a2 - D_gs4.*tan(alpha);

%Z_L_gs4= (0.25.*(f_gs4*2*pi).^2.*rho./pi./c+1j.*0.61.*rho.*f_gs4*2*pi/pi./a1).*((1+cos(alpha))/2);

Z_hole1_gs4 = 1./(-1j* ((pi*a1^2)/(rho*c))*cot(k_gs4.*(0.85.*a1+r_max_fs4)));

Z_hole2_gs4 = 1./(-1j* ((pi*a1^2)/(rho*c))*cot(k_gs4.*(0.85.*a1+r_gs4)));

Z_in_gs4 = 1./((1./ZIN1(r_max_fs4, a1, L_corr-x_fs4, 0, k_gs4, rho, c))+(1./Z_hole1_gs4));

Z_middle_gs4 = ZIN1(r_gs4, r_max_fs4, x_fs4 - D_gs4, Z_in_gs4, k_gs4, rho, c);

Z_middle2_gs4 = 1./((1./Z_middle_gs4)+(1./Z_hole2_gs4));

Z_gs4 = ZIN1(a2, r_gs4, D_gs4, Z_middle2_gs4, k_gs4, rho, c);


max_index_gs4 = find(islocalmin(db(Z_gs4)), 1);

x_gs4 = D_gs4(max_index_gs4);
r_max_gs4 = r_gs4(max_index_gs4);


figure
plot(D_gs4,db(Z_gs4));
xline(x_gs4);
xlabel("x[m]");
xlim([D_gs4(1), D_gs4(10000)]);
ylabel("Z dB");

%% PLOT OF THE RECORDER

r = linspace(a2,a1,10000);
r2 = linspace(a2,a1,10000);
x = linspace(0,L,10000);

wh = floor(a1/(x(2)-x(1)));

h1 = zeros(1,10000);
[h1_value,h1_index] = min(abs(x_fs4-x));
r(h1_index-wh : h1_index+wh) = nan;

h2 = zeros(1,10000);
[h2_value,h2_index] = min(abs(x_gs4-x));
r(h2_index-wh : h2_index+wh) = nan;

figure
plot(x,r,'k', LineWidth=2);
hold on
plot(x,-r2,'k', LineWidth=2);
yline(0, '--');
ylim([-(r(1)+0.2), r(1)+0.2]);
xline(x(h1_index),'--', 'Foot Hole');
xline(x(h2_index),'--', 'Mouth Hole');
xlabel('x[m]');
ylabel('y[m]');

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

%% QUESTION: WHY AM I DOING THIS (approximating to a normal cylindrical pipe)

n = 100000;
r_pipe = linspace(0.0001,0.2,n);

Z_L_pipe= (0.25.*(f_e4*2*pi).^2.*rho./pi./c+1j.*0.61.*rho.*f_e4*2*pi/pi./r_pipe).*((1+cos(alpha))/2);
Z_mouth_pipe = 1j* (L_c*rho./(pi.*r_pipe.^2)) .* 2.*pi.*f_e4;
Z0 = rho*c./(pi.*r_pipe.^2);
Z_pipe =  Z_mouth_pipe + Z0.* ((Z_L_pipe.*cos(k_e4.*L)+1j.*Z0.*sin(k_e4.*L))./(1j.*Z_L_pipe.*sin(k_e4.*L)+Z0.*cos(k_e4.*L)));

max_index_pipe = find(islocalmin(db(Z_pipe)), 1, 'last');

r_max_pipe = r_pipe(max_index_pipe);

figure
plot(r_pipe,db(Z_pipe));
xline(r_max_pipe);
xlabel("x[m]");
ylabel("Z dB");

