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
    a = a_vec(i)/10; Sinit = Sinit_vec(1);
    [PP1, PP2] = PP(Sinit,delta1,a);
    v1 = PP1./PP2;
    v2 = PP2;
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
    a = a_vec(i)/10; Sinit = Sinit_vec(2);
    [PP1, PP2] = PP(Sinit,delta1,a);
    v3 = PP1./PP2;
    v4 = PP2;
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

save PP.mat out_r50 out_neff

function [PP1, PP2] = PP(Sinit,delta1,a)
l = 50;
len = 100;
delta2 = 10*delta1;
Q = Sinit*delta1;
d = 1e-1;
k1 = 0.00999;
k2 = 0.999;
PP1 = zeros(len,1);
PP2 = zeros(len,1);
Mvec = logspace(-2,4,len);

for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l+1) Minit 0.1];
    [t,y] = ode15s(@EquationsPP,0:1e5:1e6,y0,odeset('RelTol',1e-5,'AbsTol',1e-10));
    PP1(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PP2(series) = sum(y(end,1:3*l+2));
end

    function ydot = EquationsPP(t,y)
        S = y(1:l+1);
        MS = y(l+2:2*l+2);
        DS = y(2*l+3:3*l+2);
        M = y(3*l+3);
        D = y(3*l+4);
        ydot = zeros(3*l+4,1);
        firstrunsumS = sum(S(1:l));
        secondrunsumS = sum(S(2:l+1));
        runsumMS = sum(MS(5:l));
        runsumDS = sum(DS(4:l));
        %% Equations for S(1), S(2), S(3), S(4), S(l+1)
        ydot(1) = Q+d*MS(1)-(a*M+delta1)*S(1)+k2*DS(1);
        ydot(2) = d*(MS(2)+DS(1))-(a*M+a*D+delta1)*S(2);
        ydot(3) = d*(MS(3)+DS(2))-(a*M+a*D+delta1)*S(3);
        ydot(4) = d*(MS(4)+DS(3))-(a*M+a*D+delta1)*S(4);
        ydot(l+1) = d*DS(l)-(a*D+delta2)*S(l+1);
        %% Equations for MS(1), MS(2), MS(3), MS(4), MS(l+1)
        ydot(l+2) = a*M*S(1)-(d+k1+delta1)*MS(1);
        ydot(l+3) = k1*MS(1)+a*M*S(2)-(d+k2+delta1)*MS(2);
        ydot(l+4) = k2*MS(2)+a*M*S(3)-(d+k2+delta1)*MS(3);
        ydot(l+5) = k2*MS(3)+a*M*S(4)-(d+k2+delta1)*MS(4);
        ydot(2*l+2) = k2*MS(l)-delta2*MS(l+1);
        %% Equations for DS(1), DS(2), DS(3), DS(l)
        ydot(2*l+3) = a*D*S(2)-(d+k2+delta1)*DS(1);
        ydot(2*l+4) = k2*DS(3)+a*D*S(3)-(d+k2+delta1)*DS(2);
        ydot(2*l+5) = k2*DS(4)+a*D*S(4)-(d+k2+delta1)*DS(3);
        ydot(3*l+2) = a*D*S(l+1)-(d+k2+delta2)*DS(l);
        %% Equations for M and D
        ydot(3*l+3) = (d+delta1)*(MS(1)+MS(2)+MS(3)+MS(4))+...
            (d+delta2)*runsumMS-a*M*firstrunsumS+delta2*MS(l+1);
        ydot(3*l+4) = (d+delta1)*(DS(1)+DS(2)+DS(3))+...
            (d+delta2)*runsumDS-a*D*secondrunsumS+k2*(DS(1)+DS(2));
        %% Equations for S(5) to S(l)
        ydot(5:l) = d*(MS(5:l)+DS(4:l-1))-(a*M+a*D+delta2)*S(5:l);
        %% Equations for MS(5) to MS(l)
        ydot(l+6:2*l+1) = k2*MS(4:l-1)+a*M*S(5:l)-(d+k2+delta2)*MS(5:l);
        %% Equations for DS(4) to DS(l-1)
        ydot(2*l+6:3*l+1) = k2*DS(5:l)+a*D*S(5:l)-(d+k2+delta2)*DS(4:l-1);
    end
end