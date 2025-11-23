global Q delta1 delta2 Moff Mon Dcat l Doff Don Mcat


Sinit = 10000;

delta1 = 2e-5;
Mon = 1e-3;
Don = 1e-4;
l = 50;
len = 10;
delta2 = 10*delta1;
Q = Sinit*delta1;
Moff = 1e-3; Doff = 1e-1;
Mcat = 0.999; Dcat = 0.999;
Dinit = 1;
DP11 = zeros(len,1);
DP12 = zeros(len,1);
Mvec = logspace(-2,4,len);
for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l) Minit Dinit];
    [t,y] = ode15s(@EquationsDP,[0 1e6],y0,odeset('RelTol',1e-6,'AbsTol',1e-8));
    DP11(series) = sum(y(end,[5:l+1,l+6:2*l+1,2*l+5:3*l+1]));
    DP12(series) = sum(y(end,1:3*l+1));
end
figure(1)
semilogx(Mvec,DP11./DP12)
figure(2)
semilogx(Mvec,DP12)