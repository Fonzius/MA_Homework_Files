clc;
clear all;
close all;

%% Ex 1.a
V = 0.1; %[m^3]
l = 0.1; %[m]
S = 100; %[m^2]
c = 343; %[m/s]
rho = 1.2; %[Kg/m^3]

r = sqrt(S/pi);

l = l + (0.61 + 0.85)*r;

Fs = 50000;

R = rho*c/S;
L = rho*l/S;
C = V/(rho*c^2);

%% Setup all the system R,L,C
subsystems = find_system('HL2','SearchDepth',2,'BlockType','SubSystem');

for ii = 1:length(subsystems)
    if contains(subsystems(ii), 'Resonator')
        set_param(char(append(subsystems(ii),'/R')), 'R', num2str(R));
        set_param(char(append(subsystems(ii),'/C')), 'C', num2str(C));
        set_param(char(append(subsystems(ii),'/L')), 'L', num2str(L));
    end
end

set_param('HL2', 'PreLoadFcn', num2str(Fs));

% Run the simulation before going on with the next parts of the code

%% Ex 1.b
f = 0:Fs/length(out.in1x1.Data):Fs-(1/length(out.in1x1.Data));
H = (abs(fft(out.out1x1.Data)./(fft(out.in1x1.Data))));

f_0 = c/(2*pi) * sqrt(S/(V*l));

maxh = 0;
maxi = 0;
for ii = 1:length(H)
    if H(ii)>maxh
        maxh = H(ii);
        maxi = ii;
    end
end

figure()
plot(f,db(H),'LineWidth',2);
xlim([0,Fs/20]);
ylim([-100,0]);
xline(f_0, 'b', 'LineWidth',1); %Theoretical resonance
xline(f(maxi), '--r', 'LineWidth',0.5); %Graph resonance
ylabel("Admittance [S]");
xlabel("f [Hz]");
legend("FRF", "Theoretical resonant Frequency", "Graph resonant Frequency");
%title("Frequency Response of the Electric-Equivalent of an Helmholtz Resonator")

error = (f(maxi)-f_0)/f_0;
f_0_graph = f(maxi);

%% Ex 2.a

H_2x2 = ((fft(out.out2x2.Data)./(fft(out.in2x2.Data))));

figure()
plot(f,db(abs(H_2x2)), 'LineWidth',2);
xlim([0,Fs/20]);
ylim([-150,0]);
ylabel("Admittance [S]");
xlabel("f [Hz]");
%title("2x2 Tree");

%% Ex 2.b

H_2x3 = (abs(fft(out.out2x3.Data)./(fft(out.in2x3.Data))));
H_3x2_3 = (abs(fft(out.out3x2_3.Data)./(fft(out.in3x2.Data))));
H_3x2i_5 = abs(fft(out.out3x2i_5.Data)./(fft(out.in3x2i.Data)));

figure()
plot(f,db(H_2x3), LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f,db(H_3x2_3), LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f,db(abs(H_2x2)), 'LineWidth',2);
xlim([0,Fs/20]);
ylim([-150,0]);
xlabel("Frequency [Hz]");
ylabel("Admittance [S]");
legend("2x3", "3x2", "2x2");
%%
figure()
plot(f,db(abs(H_3x2i_5)), 'LineWidth',2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f,db(abs(H_3x2_3)), 'LineWidth',2);
xlim([0,Fs/20]);
ylim([-150,0]);
xlabel("Frequency [Hz]");
ylabel("Admittance [S]");
legend("3x2 incomplete", "3x2");
%% Ex 2.c
%3x2 complete

H_3x2_0 = (abs(fft(out.out3x2_0.Data)./(fft(out.in3x2.Data))));
H_3x2_1 = (abs(fft(out.out3x2_1.Data)./(fft(out.in3x2.Data))));
H_3x2_2 = (abs(fft(out.out3x2_2.Data)./(fft(out.in3x2.Data))));
H_3x2_3 = (abs(fft(out.out3x2_3.Data)./(fft(out.in3x2.Data))));

figure()
plot(f, db(H_3x2_0),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f, db(H_3x2_1),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f, db(H_3x2_2),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f, db(H_3x2_3),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
xline(f_0, '--r', 'LineWidth',0.5);
xlabel("Frequency [Hz]");
ylabel("Admittance [S]");
legend("N = 0", "N = 1", "N = 2", "N = 3", "Natural frequency of single Resonator");
%title("3x2 Resonator tree response at every step");

%%

H_3x2i_0 = ((fft(out.out3x2i_0.Data)./(fft(out.in3x2i.Data))));
H_3x2i_1 = ((fft(out.out3x2i_1.Data)./(fft(out.in3x2i.Data))));
H_3x2i_2 = ((fft(out.out3x2i_2.Data)./(fft(out.in3x2i.Data))));
H_3x2i_3 = ((fft(out.out3x2i_3.Data)./(fft(out.in3x2i.Data))));
H_3x2i_4 = ((fft(out.out3x2i_4.Data)./(fft(out.in3x2i.Data))));
H_3x2i_5 = ((fft(out.out3x2i_5.Data)./(fft(out.in3x2i.Data))));

figure()
subplot(2,1,1)
plot(f, db(abs(H_3x2i_0)),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
% hold on
% plot(f, db(abs(H_3x2i_1)),LineWidth=2);
% xlim([0,Fs/20]);
% ylim([-150,0]);
% hold on
% plot(f, db(abs(H_3x2i_2)),LineWidth=2);
% xlim([0,Fs/20]);
% ylim([-150,0]);
% hold on
% plot(f, db(abs(H_3x2i_3)),LineWidth=2);
% xlim([0,Fs/20]);
% ylim([-150,0]);
% hold on
% plot(f, db(abs(H_3x2i_4)),LineWidth=2);
% xlim([0,Fs/20]);
% ylim([-150,0]);
% hold on
% plot(f, db(abs(H_3x2i_5)),LineWidth=2);
% xlim([0,Fs/20]);
% ylim([-150,0]);
% xline(f_0);
xlabel("Frequency [Hz]");
ylabel("Admittance [S]");
legend("N = 0", "N = 1", "N = 2", "N = 3", "N = 4", "N = 5");
%title("3x2 Incomplete Resonator tree response at every leaf");

subplot(2,1,2)
plot(f, unwrap(angle(H_3x2i_0)),LineWidth=2);
xlim([0,Fs/20]);
% hold on
% plot(f, unwrap(angle(H_3x2i_1)),LineWidth=2);
% xlim([0,Fs/20]);
% hold on
% plot(f, unwrap(angle(H_3x2i_2)),LineWidth=2);
% xlim([0,Fs/20]);
% hold on
% plot(f, unwrap(angle(H_3x2i_3)),LineWidth=2);
% xlim([0,Fs/20]);
% hold on
% plot(f, unwrap(angle(H_3x2i_4)),LineWidth=2);
% xlim([0,Fs/20]);
% hold on
% plot(f, unwrap(angle(H_3x2i_5)),LineWidth=2);
% xlim([0,Fs/20]);
xlabel("Frequency [Hz]");
ylabel("Phase [rad]");
legend("N = 0", "N = 1", "N = 2", "N = 3", "N = 4", "N = 5");

%%

H_test_1 = ((fft(out.out_test1.Data)./(fft(out.in_test.Data))));
H_test_2 = ((fft(out.out_test2.Data)./(fft(out.in_test.Data))));
H_test_3 = ((fft(out.out_test3.Data)./(fft(out.in_test.Data))));
H_test_4 = ((fft(out.out_test4.Data)./(fft(out.in_test.Data))));
H_test_5 = ((fft(out.out_test5.Data)./(fft(out.in_test.Data))));
H_test_6 = ((fft(out.out_test6.Data)./(fft(out.in_test.Data))));



figure()
subplot(2,1,1)
plot(f, db(abs(H_test_1)),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f, db(abs(H_test_2)),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f, db(abs(H_test_3)),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f, db(abs(H_test_4)),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f, db(abs(H_test_5)),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
hold on
plot(f, db(abs(H_test_6)),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
xline(f_0);
xlabel("Frequency [Hz]");
ylabel("Admittance [S]");
legend("N = 0", "N = 1", "N = 2", "N = 3", "N = 4", "N = 5");
%title("Incomplete Resonator tree response at every leaf");

subplot(2,1,2)
plot(f, unwrap(angle(H_test_1)),LineWidth=2);
xlim([0,Fs/20]);
hold on
plot(f, unwrap(angle(H_test_2)),LineWidth=2);
xlim([0,Fs/20]);
hold on
plot(f, unwrap(angle(H_test_3)),LineWidth=2);
xlim([0,Fs/20]);
hold on
plot(f, unwrap(angle(H_test_4)),LineWidth=2);
xlim([0,Fs/20]);
hold on
plot(f, unwrap(angle(H_test_5)),LineWidth=2);
xlim([0,Fs/20]);
hold on
plot(f, unwrap(angle(H_test_6)),LineWidth=2);
xlim([0,Fs/20]);
xlabel("Frequency [Hz]");
ylabel("Phase [rad]");
legend("N = 0", "N = 1", "N = 2", "N = 3", "N = 4", "N = 5");

%%
H_test_proof = ((fft(out.out_test_proof.Data)./(fft(out.in_test_proof.Data))));
figure()
plot(f, db(abs(H_test_proof)),LineWidth=2);
xlim([0,Fs/20]);
ylim([-150,0]);
xlabel("Frequency [Hz]");
ylabel("Admittance [S]");