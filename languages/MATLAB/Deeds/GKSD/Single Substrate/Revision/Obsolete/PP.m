function [PP1, PP2] = PP(Sinit,delta1,a)
global Q Sinit delta1 delta2 d a k1 k2 l
l = 50;
len = 100;
delta2 = 10*delta1;
Q = Sinit*delta1;
d = 1e-1;
k1 = 0.00999;
k2 = 0.999;
PP1 = zeros(len,1);
PP2 = zeros(len,1);
Mvec = logspace(-2,4,len);

for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l+1) Minit 1];
    [t,y] = ode15s(@EquationsPP,[0 1e6],y0,odeset('RelTol',1e-6,'AbsTol',1e-8));
    PP1(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PP2(series) = sum(y(end,1:3*l+2));
end
end