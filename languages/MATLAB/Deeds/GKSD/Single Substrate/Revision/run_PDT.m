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
    Mon = a_vec(i)/10; Don = a_vec(i); Sinit = Sinit_vec(1);
    [PDT1, PDT2] = PDT(Sinit,delta1,Mon,Don);
    v1 = PDT1./PDT2;
    v2 = PDT2;
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
    Mon = a_vec(i)/10; Don = a_vec(i); Sinit = Sinit_vec(2);
    [PDT1, PDT2] = PDT(Sinit,delta1,Mon,Don);
    v3 = PDT1./PDT2;
    v4 = PDT2;
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

save PDT.mat out_r50 out_neff

function [PDT1, PDT2] = PDT(Sinit,delta1,Mon,Don)
l = 50;
len = 100;
delta2 = 10*delta1;
Q = Sinit*delta1;
Moff = 1e-1;
Doff = 1e-3;
Mcat1 = 0.00999;
Mcat2 = 0.999;
Dcat = 0.999;
PDT1 = zeros(len,1);
PDT2 = zeros(len,1);
Mvec = logspace(-2,4,len);

for series = 1:len
    Minit = Mvec(series);
    y0 = [Sinit zeros(1,3*l+1) Minit 0.1];
    [t,y] = ode15s(@EquationsPDT,0:1e5:1e6,y0,odeset('RelTol',1e-5,'AbsTol',1e-10));
    PDT1(series) = sum(y(end,[5:l+1,l+6:2*l+2,2*l+6:3*l+2]));
    PDT2(series) = sum(y(end,1:3*l+2));
end

    function ydot = EquationsPDT(t,y)
        S = y(1:l+1);
        MS = y(l+2:2*l+2);
        DS = y(2*l+3:3*l+2);
        M = y(3*l+3);
        D = y(3*l+4);
        ydot = zeros(3*l+4,1);
        firstrunsumS = sum(S(1:l));
        secondrunsumS = sum(S(2:l+1));
        firstrunsumDS = sum(DS(1:l));
        secondrunsumDS = sum(DS(4:l));
        runsumMS = sum(MS(5:l));
        %% Equations for S(1), S(2), S(3), S(4), S(l+1)
        ydot(1) = Q+Moff*MS(1)-(Mon*M+delta1)*S(1)+Dcat*firstrunsumDS;
        ydot(2) = Moff*MS(2)+Doff*DS(1)-(Mon*M+Don*D+delta1)*S(2);
        ydot(3) = Moff*MS(3)+Doff*DS(2)-(Mon*M+Don*D+delta1)*S(3);
        ydot(4) = Moff*MS(4)+Doff*DS(3)-(Mon*M+Don*D+delta1)*S(4);
        ydot(l+1) = Doff*DS(l)-(Don*D+delta2)*S(l+1);
        %% Equations for MS(1), MS(2), MS(3), MS(4), MS(l+1)
        ydot(l+2) = Mon*M*S(1)-(Moff+Mcat1+delta1)*MS(1);
        ydot(l+3) = Mcat1*MS(1)+Mon*M*S(2)-(Moff+Mcat2+delta1)*MS(2);
        ydot(l+4) = Mcat2*MS(2)+Mon*M*S(3)-(Moff+Mcat2+delta1)*MS(3);
        ydot(l+5) = Mcat2*MS(3)+Mon*M*S(4)-(Moff+Mcat2+delta1)*MS(4);
        ydot(2*l+2) = Mcat2*MS(l)-delta2*MS(l+1);
        %% Equations for DS(1), DS(2), DS(3), DS(4), DS(l+1)
        ydot(2*l+3) = Don*D*S(2)-(Doff+Dcat+delta1)*DS(1);
        ydot(2*l+4) = Don*D*S(3)-(Doff+Dcat+delta1)*DS(2);
        ydot(2*l+5) = Don*D*S(4)-(Doff+Dcat+delta1)*DS(3);
        ydot(2*l+6) = Don*D*S(5)-(Doff+Dcat+delta2)*DS(4);
        ydot(3*l+2) = Don*D*S(l+1)-(Doff+Dcat+delta2)*DS(l);
        %% Equations for M and D
        ydot(3*l+3) = (Moff+delta1)*(MS(1)+MS(2)+MS(3)+MS(4))+...
            (Moff+delta2)*runsumMS-Mon*M*firstrunsumS+delta2*MS(l+1);
        ydot(3*l+4) = (Doff+Dcat+delta1)*(DS(1)+DS(2)+DS(3))+...
            (Doff+Dcat+delta2)*secondrunsumDS-Don*D*secondrunsumS;
        %% Equations for S(5) to S(l)
        ydot(5:l) = Moff*MS(5:l)+Doff*DS(4:l-1)-(Mon*M+Don*D+delta2)*S(5:l);
        %% Equations for MS(5) to MS(l)
        ydot(l+6:2*l+1) = Mcat2*MS(4:l-1)+Mon*M*S(5:l)-(Moff+Mcat2+delta2)*MS(5:l);
        %% Equations for DS(5) to DS(l)
        ydot(2*l+7:3*l+1) = Don*D*S(6:l)-(Doff+Dcat+delta2)*DS(5:l-1);
    end
end