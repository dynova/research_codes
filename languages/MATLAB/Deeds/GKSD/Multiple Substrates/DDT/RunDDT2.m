global Q1 Q2 r1 r2 d a k l
l = 50;
len = 100;
S1init = 100;
r1 = 2e-5;
r2 = 10*r1;
Q1 = S1init*r1;
d = 1e-3;
a = 1e-3;
k = 0.999;
DDT1 = zeros(len,1);
DDT2 = zeros(len,1);
S2initvec = logspace(1,5,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    S2init = S2initvec(series);
    Q2 = S2init*r1;
    y0 = [S1init zeros(1,3*l) S2init zeros(1,3*l) 1 1];
    [t,y] = ode15s(@EquationsDDT2,0:1e5:1e6,y0,options);
    DDT1(series) = sum(y(end,[5:l+1,l+6:2*l+1,2*l+5:3*l+1]));
    DDT2(series) = sum(y(end,1:3*l+1));
end
save DDT DDT1 DDT2