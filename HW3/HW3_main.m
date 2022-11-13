% HW3_main

n_L=1:100; %length interval [m]
e1=zeros(length(n_L),1);
delta=zeros(length(n_L),1);

for i=1:length(n_L)
    [e1(i,1),delta(i,1)]=ERROR_e(n_L(i));
end
plot(delta,e1)



