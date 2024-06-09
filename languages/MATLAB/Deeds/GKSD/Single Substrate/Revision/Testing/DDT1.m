global Q delta1 delta2 aM aD dM dD kcatM kcatD l


Sinit = 100;

delta1 = 2e-5;
aM = 1e-3/10;
aD = 1e-3/10;
l = 50;
len = 100;
delta2 = 10*delta1;
Q = Sinit*delta1;
dM = 1e-3; dD = 1e-3;
kcatM = 0.999; kcatD = 0.999;
Dinit = 1;
DDT11 = zeros(len,1);
DDT12 = zeros(len,1);
Mvec = logspace(-2,4,len);
for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l) Minit Dinit];
    [t,y] = ode15s(@EquationsDDT,[0 1e6],y0,odeset('RelTol',1e-6,'AbsTol',1e-8));
    DDT11(series) = sum(y(end,[5:l+1,l+6:2*l+1,2*l+5:3*l+1]));
    DDT12(series) = sum(y(end,1:3*l+1));
end
figure(1)
semilogx(Mvec,DDT11./DDT12)
figure(2)
semilogx(Mvec,DDT12)