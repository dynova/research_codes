function [S1mod, S1tot, S2mod, S2tot] = Task1Prog5
global Q1 Q2 delta1 delta2 Moff Mon Mcat Doff Don Dcat

E3vec = logspace(-2,2,100);
lenE3 = length(E3vec);
S11plot = zeros(lenE3,1);
S12plot = zeros(lenE3,1);
MS1plot = zeros(lenE3,1);
DS1plot = zeros(lenE3,1);
S21plot = zeros(lenE3,1);
S22plot = zeros(lenE3,1);
MS2plot = zeros(lenE3,1);
DS2plot = zeros(lenE3,1);
delta1 = 2e-5;
delta2 = 2e-4;
Mon = 1e-3;
Moff = 0.001;
Mcat = 0.999;
Don = 1e-3;
Doff = 0.001;
Dcat = 0.999;
Q1 = 1;
Q2 = 1e-0;
tf = 1e8;
time = 0:tf/1E3:tf;

for E3loop = 1:lenE3
    E30 = E3vec(E3loop);
    y0 = [Q1/delta1 0 0 0 Q2/delta1 0 0 0 E30 1];
    [t,y] = ode15s(@Equations,time,y0,odeset('RelTol',1e-8,'AbsTol',1e-8));
    S11plot(E3loop,1) = y(lenE3,1);
    S12plot(E3loop,1) = y(lenE3,2);
    MS1plot(E3loop,1) = y(lenE3,3);
    DS1plot(E3loop,1) = y(lenE3,4);
    S21plot(E3loop,1) = y(lenE3,5);
    S22plot(E3loop,1) = y(lenE3,6);
    MS2plot(E3loop,1) = y(lenE3,7);
    DS2plot(E3loop,1) = y(lenE3,8);
    S1mod = S12plot + MS1plot + DS1plot;
    S1tot = S11plot + S12plot + MS1plot + DS1plot;
    S2mod = S22plot + MS2plot + DS2plot;
    S2tot = S21plot + S22plot + MS2plot + DS2plot;
end
    function ydot = Equations(t,y)
        S11 = y(1);
        S12 = y(2);
        MS1 = y(3);
        DS1 = y(4);
        S21 = y(5);
        S22 = y(6);
        MS2 = y(7);
        DS2 = y(8);
        M = y(9);
        D = y(10);
        ydot = zeros(10,1);
        ydot(1) = Q1+Moff*MS1-(Mon*M+delta1)*S11+Dcat*DS1;
        ydot(2) = Doff*DS1-(Don*D+delta2)*S12+Mcat*MS1;
        ydot(3) = Mon*M*S11-(Moff+Mcat+delta1)*MS1;
        ydot(4) = Don*D*S12-(Doff+Dcat+delta2)*DS1;
        ydot(5) = Q2+Moff*MS2-(Mon*M+delta1)*S21+Dcat*DS2;
        ydot(6) = Doff*DS2-(Don*D+delta2)*S22+Mcat*MS2;
        ydot(7) = Mon*M*S21-(Moff+Mcat+delta1)*MS2;
        ydot(8) = Don*D*S22-(Doff+Dcat+delta2)*DS2;
        ydot(9) = (Moff+Mcat+delta1)*(MS1+MS2)-Mon*M*(S11+S21);
        ydot(10) = (Doff+Dcat+delta2)*(DS1+DS2)-Don*D*(S12+S22);
    end
end