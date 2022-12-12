clc
clear
close

peaknumbers=10;

%% Change the length of exponential horn
ll=0.15:0.1:1.05; %ori=0.35
X_compound=zeros(length(ll),peaknumbers);
F_compound=zeros(length(ll),peaknumbers);
for ii=1:length(ll)

[x_compound, f_compound]= ChangeLengthhorn(ll(ii),peaknumbers);
X_compound(ii,:)=x_compound;
F_compound(ii,:)=f_compound;
end

Strlegend='';
figure()
for jj=1:length(ll)
plot(X_compound(jj,:),F_compound(jj,:),'--',LineWidth=1.5)
nlegend(jj)=strcat(string(ll(jj)),'[m]');
% Strlegend=strcat(nlegend);

% legend(jj)
hold on
end

hold on

xx=1:peaknumbers;
yy=1:peaknumbers;
yy2=1:2:2*peaknumbers;
plot(xx,yy,'r',xx,yy2,'m',LineWidth=2)

legend(nlegend(1),nlegend(2),nlegend(3),nlegend(4),nlegend(5),...
    nlegend(6),nlegend(7),nlegend(8),nlegend(9),nlegend(10),... 
    'harmonic','harmonic-odd');
xlim([0 14])
xlabel('Peaknumber')
ylabel('Mulipules')

%% Change the length of the cylinder
ll=0.001:0.3:2.25; %ori=0.85
X_compound=zeros(length(ll),peaknumbers);
F_compound=zeros(length(ll),peaknumbers);
for ii=1:length(ll)

[x_compound, f_compound]= ChangeLengthcylinder(ll(ii),peaknumbers);
X_compound(ii,:)=x_compound;
F_compound(ii,:)=f_compound;
end

Strlegend='';
figure()
for jj=1:length(ll)
plot(X_compound(jj,:),F_compound(jj,:),'--',LineWidth=1.5)
nlegend(jj)=strcat(string(ll(jj)),'[m]');
% Strlegend=strcat(nlegend);

% legend(jj)
hold on
end

hold on

xx=1:peaknumbers;
yy=1:peaknumbers;
yy2=1:2:2*peaknumbers;
plot(xx,yy,'r',xx,yy2,'m',LineWidth=2)
legend()
xlabel('Peaknumber')
ylabel('Mulipules')
legend(nlegend(1),nlegend(2),nlegend(3),nlegend(4),nlegend(5),...
    nlegend(6),nlegend(7),nlegend(8),... 
    'harmonic','harmonic-odd');
xlim([0 14])
