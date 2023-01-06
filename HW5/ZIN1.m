function Zin = ZIN1 (r1,r2,L,Zl,k,rho,c)

S1=pi.*r1.^2;
S2=pi.*r2.^2;
x1=r1.*L/(r2-r1);
x2=x1+L;
theta1=atan(k.*x1);
theta2=atan(k.*x2);

Z_n=1j.*Zl.*(sin(k.*L-theta2)./sin(theta2))+rho.*c./S2.*sin(k.*L);
Z_d=Zl.*(sin(k.*L+theta1-theta2)./(sin(theta1).*sin(theta2)))-...
    (1j.*rho.*c./S2).*(sin(k.*L+theta1)./sin(theta1));

Zin=rho.*c./S1.*Z_n./Z_d;


