function ydot = EquationsPDS(t,y)
global Q delta1 delta2 Doff Don Dcat Mcat1 Mcat2 l Mon Moff
S = y(1:l+1);
MS = y(l+2:2*l+2);
DS = y(2*l+3:3*l+2);
M = y(3*l+3);
D = y(3*l+4);
ydot = zeros(3*l+4,1);
firstrunsumS = sum(S(1:l));
secondrunsumS = sum(S(2:l+1));
runsumMS = sum(MS(5:l));
runsumDS = sum(DS(5:l));
%% Equations for S(1), S(2), S(3), S(4), S(l+1)
ydot(1) = Q+Moff*MS(1)-(Mon*M+delta1)*S(1)+Dcat*DS(1);
ydot(2) = Moff*MS(2)+Doff*DS(1)-(Mon*M+Don*D+delta1)*S(2)+Dcat*DS(2);
ydot(3) = Moff*MS(3)+Doff*DS(2)-(Mon*M+Don*D+delta1)*S(3)+Dcat*DS(3);
ydot(4) = Moff*MS(4)+Doff*DS(3)-(Mon*M+Don*D+delta1)*S(4)+Dcat*DS(4);
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
ydot(2*l+6) = Don*D*S(5)-(Doff+Dcat+delta1)*DS(4);
ydot(3*l+2) = Don*D*S(l+1)-(Doff+Dcat+delta2)*DS(l);
%% Equations for M and D
ydot(3*l+3) = (Moff+delta1)*(MS(1)+MS(2)+MS(3)+MS(4))+...
    (Moff+delta2)*runsumMS-Mon*M*firstrunsumS+delta2*MS(l+1);
ydot(3*l+4) = (Doff+Dcat+delta1)*(DS(1)+DS(2)+DS(3)+DS(4))+...
    (Doff+Dcat+delta2)*runsumDS-Don*D*secondrunsumS;
%% Equations for S(5) to S(l)
ydot(5:l) = Moff*MS(5:l)+Doff*DS(4:l-1)-(Mon*M+Don*D+delta2)*S(5:l)+Dcat*DS(5:l);
%% Equations for MS(5) to MS(l)
ydot(l+6:2*l+1) = Mcat2*MS(4:l-1)+Mon*M*S(5:l)-(Moff+Mcat2+delta2)*MS(5:l);
%% Equations for DS(5) to DS(l)
ydot(2*l+7:3*l+1) = Don*D*S(6:l)-(Doff+Dcat+delta2)*DS(5:l-1);
end