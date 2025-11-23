rng shuffle;
Sinit_vec = [1e2 1e4];
d1_base = 2e-5;
d1_left = log10(d1_base) - 1;
d1_right = log10(d1_base) + 1;
a_base = 1e-3;
a_left = log10(a_base) - 1;
a_right = log10(a_base) + 1;
n = 100;
d1_vec = 10.^(d1_left + (d1_right-d1_left) * rand(1,n));
a_vec = 10.^(a_left + (a_right-a_left) * rand(1,n));
count = 0;
[out1, out2, out3, out4] = deal(zeros(n,2));

for i = 1:n
    delta1 = d1_vec(i);
    aM = a_vec(i); aD = a_vec(i); Sinit = Sinit_vec(1);
    [DDS1, DDS2] = DDS(Sinit,delta1,aM,aD);
    v1 = DDS1./DDS2;
    v2 = DDS2;
    [~,v1_r10] = min(abs(v1-(0.9*min(v1)+0.1*max(v1))));
    [~,v2_r10] = min(abs(v2-(0.9*min(v2)+0.1*max(v2))));
    [~,v1_r50] = min(abs(v1-(0.5*min(v1)+0.5*max(v1))));
    [~,v2_r50] = min(abs(v2-(0.5*min(v2)+0.5*max(v2))));
    [~,v1_r90] = min(abs(v1-(0.1*min(v1)+0.9*max(v1))));
    [~,v2_r90] = min(abs(v2-(0.1*min(v2)+0.9*max(v2))));
    v1_neff = log(81)/log(v1_r90/v1_r10);
    v2_neff = log(81)/log(v2_r10/v2_r90);
    out1(i,:) = [v2_r50, v1_r50];
    out2(i,:) = [v2_neff, v1_neff];
    count = count + 1;
    disp(count)
end
for i = 1:n
    delta1 = d1_vec(i);
    aM = a_vec(i); aD = a_vec(i); Sinit = Sinit_vec(2);
    [DDS1, DDS2] = DDS(Sinit,delta1,aM,aD);
    v3 = DDS1./DDS2;
    v4 = DDS2;
    [~,v3_r10] = min(abs(v3-(0.9*min(v3)+0.1*max(v3))));
    [~,v4_r10] = min(abs(v4-(0.9*min(v4)+0.1*max(v4))));
    [~,v3_r50] = min(abs(v3-(0.5*min(v3)+0.5*max(v3))));
    [~,v4_r50] = min(abs(v4-(0.5*min(v4)+0.5*max(v4))));
    [~,v3_r90] = min(abs(v3-(0.1*min(v3)+0.9*max(v3))));
    [~,v4_r90] = min(abs(v4-(0.1*min(v4)+0.9*max(v4))));
    v3_neff = log(81)/log(v3_r90/v3_r10);
    v4_neff = log(81)/log(v4_r10/v4_r90);
    out3(i,:) = [v4_r50, v3_r50];
    out4(i,:) = [v4_neff, v3_neff];
    count = count + 1;
    disp(count)
end
out_r50 = out3./out1;
out_neff = out4./out2;

save DDS.mat out_r50 out_neff

function [DDS1, DDS2] = DDS(Sinit,delta1,aM,aD)
l = 50;
len = 100;
delta2 = 10*delta1;
Q = Sinit*delta1;
dM = 1e-3; dD = 1e-3;
kcatM = 0.999; kcatD = 0.999;
Dinit = 0.1;
DDS1 = zeros(len,1);
DDS2 = zeros(len,1);
Mvec = logspace(-2,4,len);

for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l) Minit Dinit];
    [t,y] = ode15s(@EquationsDDS,0:1e5:1e6,y0,odeset('RelTol',1e-5,'AbsTol',1e-10));
    DDS1(series) = sum(y(end,[5:l+1,l+6:2*l+1,2*l+5:3*l+1]));
    DDS2(series) = sum(y(end,1:3*l+1));
end

    function ydot = EquationsDDS(t,y)
        S = y(1:l+1);
        MS = y(l+2:2*l+1);
        DS = y(2*l+2:3*l+1);
        M = y(3*l+2);
        D = y(3*l+3);
        ydot = zeros(3*l+3,1);
        firstrunsumS = sum(S(1:l));
        secondrunsumS = sum(S(2:l+1));
        runsumMS = sum(MS(5:l));
        runsumDS = sum(DS(4:l));
        %% Equations for S0, S1, S2, S3, S_(l)
        ydot(1) = Q-delta1*S(1)+dM*MS(1)-aM*M*S(1)+kcatD*DS(1);
        ydot(2) = -S(2)*(delta1+aM*M+aD*D)+kcatM*MS(1)+kcatD*DS(2)+dD*DS(1)+dM*MS(2);
        ydot(3) = -S(3)*(delta1+aM*M+aD*D)+kcatM*MS(2)+kcatD*DS(3)+dD*DS(2)+dM*MS(3);
        ydot(4) = -S(4)*(delta1+aM*M+aD*D)+kcatM*MS(3)+kcatD*DS(4)+dD*DS(3)+dM*MS(4);
        ydot(l+1) = kcatM*MS(l) - (aD*D+delta2)*S(l+1) + dD*DS(l);
        %% Equations for MS0, MS1, MS2, MS3
        ydot(l+2) = aM*M*S(1)-(dM+kcatM+delta1)*MS(1);
        ydot(l+3) = aM*M*S(2)-(dM+kcatM+delta1)*MS(2);
        ydot(l+4) = aM*M*S(3)-(dM+kcatM+delta1)*MS(3);
        ydot(l+5) = aM*M*S(4)-(dM+kcatM+delta1)*MS(4);
        %% Equations for DS1, DS2, DS3
        ydot(2*l+2) = aD*D*S(2)-(dD+kcatD+delta1)*DS(1);
        ydot(2*l+3) = aD*D*S(3)-(dD+kcatD+delta1)*DS(2);
        ydot(2*l+4) = aD*D*S(4)-(dD+kcatD+delta1)*DS(3);
        %% Equations for M and D
        ydot(3*l+2) = (dM+kcatM+delta1)*(MS(1)+MS(2)+MS(3)+MS(4)) + (dM+kcatM+delta2)*runsumMS - ...
            aM*M*firstrunsumS;
        ydot(3*l+3) = (dD+kcatD+delta1)*(DS(1)+DS(2)+DS(3)) + (dD+kcatD+delta2)*runsumDS - ...
            aD*D*secondrunsumS;
        %% Equations for S4, S5, ..., S_(l-1) and MS4, MS5, ..., MS_(l-1)
        ydot(5:l) = -S(5:l)*(delta2+aM*M+aD*D)+kcatM*MS(4:l-1)+kcatD*DS(5:l)+ ...
            dD*DS(4:l-1)+dM*MS(5:l);
        ydot(l+6:2*l+1)=aM*M*S(5:l)-(dM+kcatM+delta2)*MS(5:l);
        %% Equations for DS4, DS5, ..., DS_(l)
        ydot(2*l+5:3*l+1)=aD*D*S(5:l+1)-(dD+kcatD+delta2)*DS(4:l);
    end
end