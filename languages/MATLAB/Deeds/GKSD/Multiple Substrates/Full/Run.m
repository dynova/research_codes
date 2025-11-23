global Q1 Q2 delta1 delta2 Moff Mon Mcat Doff Don Dcat
len = 100;
S1init = 1e4;
delta1 = 2e-5;
delta2 = 10*delta1;
Q1 = S1init*delta1;
Moff = 1e-3;
Mon = 1e-3;
Mcat = 0.999;
Doff = 1e-3;
Don = 1e-3;
Dcat = 0.999;
A = zeros(len,1);
S2initvec = logspace(1,9,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    S2init = S2initvec(series);
    Q2 = S2init*delta1;
    y0 = [S1init zeros(1,3) S2init zeros(1,3) 0.1 0.1];
    [t,y] = ode15s(@Equations,0:1e7:1e8,y0,options);
    A(series) = sum(y(end,1:4));
end
save Full.mat A;