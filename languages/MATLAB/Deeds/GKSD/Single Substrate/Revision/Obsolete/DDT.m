function [DDT1, DDT2] = DDT(Sinit,delta1,aM,aD)

global Q delta1 delta2 dM dD aM aD kcatM kcatD l
l = 50;
len = 100;
delta2 = 10*delta1;
Q = Sinit*delta1;
dM = 1e-3; dD = 1e-3;
kcatM = 0.999; kcatD = 0.999;
Dinit = 1;
DDT1 = zeros(len,1);
DDT2 = zeros(len,1);
Mvec = logspace(-2,4,len);

for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l) Minit Dinit];
    [t,y] = ode15s(@EquationsDDT,[0 1e6],y0,odeset('RelTol',1e-6,'AbsTol',1e-8));
    DDT1(series) = sum(y(end,[5:l+1,l+6:2*l+1,2*l+5:3*l+1]));
    DDT2(series) = sum(y(end,1:3*l+1));
end
end