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
PP41 = zeros(len,1);
PP42 = zeros(len,1);
Mvec = logspace(-2,1,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l+1) Minit 0.1];
    [t,y] = ode15s(@EquationsPP,[0 1e6],y0,options);
    PP41(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PP42(series) = sum(y(end,1:3*l+2));
end
save PP4 PP41 PP42