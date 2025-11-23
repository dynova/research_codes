global Q1 Q2 delta1 delta2 Doff Don Dcat Mcat1 Mcat2 l Mon Moff
S1init = 1e4;
Mcat1 = 0.00999;
l = 50;
len = 100;
delta1 = 2e-5;
delta2 = 10*delta1;
Q1 = S1init*delta1;
Moff = 1e-1;
Mon = 1e-4;
Mcat2 = 0.999;
Doff = 1e-3;
Don = 1e-3;
Dcat = 0.999;
B = zeros(1,len);
S2initvec = logspace(1,9,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for ind = 1:len
    S2init = S2initvec(ind);
    Q2 = S2init*delta1;
    y0 = [S1init zeros(1,3*l+1) S2init zeros(1,3*l+1) 100 1];
    [t,y] = ode15s(@EquationsPDS,0:1e5:1e6,y0,options);
    B(ind) = sum(y(end,1:3*l+2));
end
save PDS.mat B;