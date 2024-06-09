global Q delta1 delta2 Moff Mon Dcat l Doff Don Mcat
l = 50;
len = 100;
Sinit = 10000;
delta1 = 2e-5;
delta2 = 10*delta1;
Q = Sinit*delta1;
Moff = 1e-3;
Doff = 1e-1;
Mon = 1e-3;
Don = 1e-4;
Mcat = 0.999;
Dcat = 0.999;
DP41 = zeros(len,1);
DP42 = zeros(len,1);
Mvec = logspace(-2,1,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l) Minit 0.1];
    [t,y] = ode15s(@EquationsDP,0:1e5:1e6,y0,options);
    DP41(series) = sum(y(end,[5:l+1,l+6:2*l+1,2*l+5:3*l+1]));
    DP42(series) = sum(y(end,1:3*l+1));
end
save DP4 DP41 DP42