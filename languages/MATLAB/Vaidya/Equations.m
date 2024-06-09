function ydot = Equations(t,y)

global Lambda mu mu_A mu_T beta tau eta_1 eta_2 lambda_T rho_1 rho_2 gamma

S = y(1);
T = y(2);
H = y(3);
C = y(4);
C1 = y(5); 
C2 = y(6);
CM1 = y(7); 
CM2 = y(8);

N = S + T + H + C + C1 + C2 + CM1 + CM2;
ydot = zeros(8,1);

ydot(1)=Lambda-mu.*S-beta.*(H+C+C1+C2).*(S./N)-tau.*(T+C).*(S./N);
ydot(2)=tau.*(T+C).*(S./N)-beta.*(H+C+C1+C2).*(T./N)-(mu+mu_T).*T;
ydot(3)=beta.*(H+C+C1+C2).*(S./N)-tau.*(T+C).*(H./N)-(mu+mu_A).*H;
ydot(4)=beta.*(H+C+C1+C2).*(T./N)+tau.*(T+C).*(H./N)-(mu+mu_A+mu_T+lambda_T).*C;
ydot(5)=lambda_T.*C-(mu+mu_A+rho_1+eta_1).*C1;
ydot(6)=rho_1.*C1-(mu+mu_A+rho_2+eta_2).*C2;
ydot(7)=eta_1.*C1-(mu+rho_1+gamma).*CM1;
ydot(8)=eta_2.*C2+(rho_1).*CM1-(mu+rho_2+gamma.*(rho_1)./(rho_1+rho_2)).*CM2;
end