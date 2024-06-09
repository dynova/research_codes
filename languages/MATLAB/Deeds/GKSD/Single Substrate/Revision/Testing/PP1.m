global Q Sinit delta1 delta2 d a k1 k2 l


Sinit = 100;

delta1 = 2e-5;
a = 1e-4;
l = 50;
len = 10;
delta2 = 10*delta1;
Q = Sinit*delta1;
d = 1e-1;
k1 = 0.00999;
k2 = 0.999;
PP11 = zeros(len,1);
PP12 = zeros(len,1);
Mvec = logspace(-2,4,len);
for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l+1) Minit 1];
    [t,y] = ode15s(@EquationsPP,[0 1e6],y0,odeset('RelTol',1e-6,'AbsTol',1e-8));
    PP11(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PP12(series) = sum(y(end,1:3*l+2));
end
figure(1)
semilogx(Mvec,PP11./PP12)
figure(2)
semilogx(Mvec,PP12)