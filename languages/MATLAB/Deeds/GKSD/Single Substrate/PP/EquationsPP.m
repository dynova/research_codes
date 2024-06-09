function ydot = EquationsPP(t,y)
global Q delta1 delta2 d a k1 k2 l
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