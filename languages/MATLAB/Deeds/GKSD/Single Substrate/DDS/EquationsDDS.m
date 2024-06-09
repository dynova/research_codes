function ydot = EquationsDDS(t,y)
global Q delta1 delta2 dM dD aM aD kcatM kcatD l
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