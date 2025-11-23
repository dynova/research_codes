function [DP1, DP2] = DP(Sinit,delta1,Mon,Don)

global Q delta1 delta2 Moff Mon Dcat l Doff Don Mcat
l = 50;
len = 100;
delta2 = 10*delta1;
Q = Sinit*delta1;
Moff = 1e-3; Doff = 1e-1;
Mcat = 0.999; Dcat = 0.999;
Dinit = 1;
DP1 = zeros(len,1);
DP2 = zeros(len,1);
Mvec = logspace(-2,4,len);

for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l) Minit Dinit];
    [t,y] = ode15s(@EquationsDP,[0 1e6],y0,odeset('RelTol',1e-6,'AbsTol',1e-8));
    DP1(series) = sum(y(end,[5:l+1,l+6:2*l+1,2*l+5:3*l+1]));
    DP2(series) = sum(y(end,1:3*l+1));
end
end