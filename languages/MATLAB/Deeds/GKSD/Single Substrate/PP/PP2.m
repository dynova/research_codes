global Q delta1 delta2 d a k1 k2 l
l = 50;
len = 100;
Sinit = 10000;
delta1 = 2e-5;
delta2 = 10*delta1;
Q = Sinit*delta1;
d = 1e-1;
a = 1e-4;
k1 = 0.00999;
k2 = 0.999;
PP21 = zeros(len,1);
PP22 = zeros(len,1);
Mvec = logspace(-2,4,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l+1) Minit 0.1];
    [t,y] = ode15s(@EquationsPP,[0 1e6],y0,options);
    PP21(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PP22(series) = sum(y(end,1:3*l+2));
end
save PP2 PP21 PP22