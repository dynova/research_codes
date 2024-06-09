global Q1 Q2 delta1 delta2 Moff Mon Dcat l Doff Don Mcat
l = 50;
len = 10;
S1init = 100;
delta1 = 2e-5;
delta2 = 10*delta1;
Q1 = S1init*delta1;
Doff = 1e-1;
Don = 1e-4;
Dcat = 0.999;
Mcat = 0.999;
Mon = 1e-3;
Moff = 1e-3;
DP1 = zeros(len,1);
DP2 = zeros(len,1);
S2initvec = logspace(1,5,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    S2init = S2initvec(series);
    Q2 = S2init*delta1;
    y0 = [S1init zeros(1,3*l) S2init zeros(1,3*l) 1 1];
    [t,y] = ode15s(@EquationsDP2,0:1e5:1e6,y0,options);
    DP1(series) = sum(y(end,[5:l+1,l+6:2*l+1,2*l+5:3*l+1]));
    DP2(series) = sum(y(end,1:3*l+1));
end
save DPM DP1 DP2