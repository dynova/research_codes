global Q1 Q2 delta1 delta2 Doff Don Dcat Mcat1 Mcat2 l Mon Moff
l = 500;
len = 100;
S1init = 100;
delta1 = 2e-5;
delta2 = 10*delta1;
Q1 = S1init*delta1;
Doff = 1e-3;
Don = 1e-3;
Dcat = 0.999;
Mcat1 = 0.00999;
Mcat2 = 0.999;
Mon = 1e-4;
Moff = 1e-1;
PDT1 = zeros(len,1);
PDT2 = zeros(len,1);
S2initvec = logspace(1,5,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    S2init = S2initvec(series);
    Q2 = S2init*delta1;
    y0 = [S1init zeros(1,3*l+1) S2init zeros(1,3*l+1) 1 1];
    [t,y] = ode15s(@EquationsPDT2,0:1e5:1e6,y0,options);
    PDT1(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PDT2(series) = sum(y(end,1:3*l+2));
end
save PDT PDT1 PDT2