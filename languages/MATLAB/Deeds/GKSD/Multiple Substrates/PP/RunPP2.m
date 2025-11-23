global Q1 Q2 delta1 delta2 d a k1 k2 l
l = 500;
len = 100;
S1init = 100;
delta1 = 2e-5;
delta2 = 10*delta1;
Q1 = S1init*delta1;
d = 1e-1;
a = 1e-4;
k1 = 0.00999;
k2 = 0.999;
PP1 = zeros(len,1);
PP2 = zeros(len,1);
S2initvec = logspace(1,5,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    S2init = S2initvec(series);
    Q2 = S2init*delta1;
    y0 = [S1init zeros(1,3*l+1) S2init zeros(1,3*l+1) 1 1];
    [t,y] = ode15s(@EquationsPP2,0:1e5:1e6,y0,options);
    PP1(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PP2(series) = sum(y(end,1:3*l+2));
end
save PP PP1 PP2