global Q delta1 delta2 Doff Don Dcat Mcat1 Mcat2 l Mon Moff
l = 50;
len = 100;
Sinit = 10000;
delta1 = 2e-5;
delta2 = 10*delta1;
Q = Sinit*delta1;
Moff = 1e-1;
Doff = 1e-3;
Mon = 1e-4;
Don = 1e-3;
Mcat1 = 0.00999;
Mcat2 = 0.999;
Dcat = 0.999;
PDT41 = zeros(len,1);
PDT42 = zeros(len,1);
Mvec = logspace(-2,1,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l+1) Minit 0.1];
    [t,y] = ode15s(@EquationsPDT,0:1e5:1e6,y0,options);
    PDT41(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PDT42(series) = sum(y(end,1:3*l+2));
end
save PDT4 PDT41 PDT42