global Q Sinit delta1 delta2 Doff Dcat Mcat1 Mcat2 l Moff Mon Don


Sinit = 10000;

delta1 = 2e-5;
Mon = 1e-4;
Don = 1e-3;
l = 50;
len = 100;
delta2 = 10*delta1;
Q = Sinit*delta1;
Moff = 1e-1;
Doff = 1e-3;
Mcat1 = 0.00999;
Mcat2 = 0.999;
Dcat = 0.999;
PDT11 = zeros(len,1);
PDT12 = zeros(len,1);
Mvec = logspace(-2,4,len);
for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l+1) Minit 1];
    [t,y] = ode15s(@EquationsPDT,[0 1e6],y0,odeset('RelTol',1e-6,'AbsTol',1e-8));
    PDT11(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PDT12(series) = sum(y(end,1:3*l+2));
end
figure(1)
semilogx(Mvec,PDT11./PDT12)
figure(2)
semilogx(Mvec,PDT12)