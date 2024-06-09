global Q delta1 delta2 dM dD aM aD kcatM kcatD l
l = 50;
len = 100;
Sinit = 100;
delta1 = 2e-5;
delta2 = 10*delta1;
Q = Sinit*delta1;
dM = 1e-3;
dD = 1e-3;
aM = 1e-3;
aD = 1e-3;
kcatM = 0.999;
kcatD = 0.999;
DDS31 = zeros(len,1);
DDS32 = zeros(len,1);
Mvec = logspace(-2,1,len);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l) Minit 0.1];
    [t,y] = ode15s(@EquationsDDS,0:1e5:1e6,y0,options);
    DDS31(series) = sum(y(end,[5:l+1,l+6:2*l+1,2*l+5:3*l+1]));
    DDS32(series) = sum(y(end,1:3*l+1));
end
save DDS3.mat DDS31 DDS32