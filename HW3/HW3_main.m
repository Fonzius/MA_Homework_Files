% HW3_main

n_L=100; %length interval [m] 1:100

%% a)
e1=zeros(length(n_L),1);
delta=zeros(length(n_L),1);

for i=1:length(n_L)
    [e1(i,1),delta(i,1)]=ERROR_e(n_L(i));
end
figure();
plot(delta,e1)

% %% b)
% 
% e2=zeros(length(n_L),1);
% delta=zeros(length(n_L),1);
% 
% for i=





