function ydot = EquationsPP2(t,y)
global Q1 Q2 delta1 delta2 d a k1 k2 l
S1 = y(1:l+1);
MS1 = y(l+2:2*l+2);
DS1 = y(2*l+3:3*l+2);
S2 = y(3*l+3:4*l+3);
MS2 = y(4*l+4:5*l+4);
DS2 = y(5*l+5:6*l+4);
M = y(6*l+5);
D = y(6*l+6);
ydot = zeros(6*l+6,1);
firstrunsumS1 = sum(S1(1:l));
secondrunsumS1 = sum(S1(2:l+1));
runsumMS1 = sum(MS1(5:l));
runsumDS1 = sum(DS1(4:l));
firstrunsumS2 = sum(S2(1:l));
secondrunsumS2 = sum(S2(2:l+1));
runsumMS2 = sum(MS2(5:l));
runsumDS2 = sum(DS2(4:l));
%% Equations for S(1), S(2), S(3), S(4), S(l+1)
ydot(1) = Q1+d*MS1(1)-(a*M+delta1)*S1(1)+k2*DS1(1);
ydot(2) = d*(MS1(2)+DS1(1))-(a*M+a*D+delta1)*S1(2);
ydot(3) = d*(MS1(3)+DS1(2))-(a*M+a*D+delta1)*S1(3);
ydot(4) = d*(MS1(4)+DS1(3))-(a*M+a*D+delta1)*S1(4);
ydot(5:l) = d*(MS1(5:l)+DS1(4:l-1))-(a*M+a*D+delta2)*S1(5:l);
ydot(l+1) = d*DS1(l)-(a*D+delta2)*S1(l+1);
ydot(l+2) = a*M*S1(1)-(d+k1+delta1)*MS1(1);
ydot(l+3) = k1*MS1(1)+a*M*S1(2)-(d+k2+delta1)*MS1(2);
ydot(l+4) = k2*MS1(2)+a*M*S1(3)-(d+k2+delta1)*MS1(3);
ydot(l+5) = k2*MS1(3)+a*M*S1(4)-(d+k2+delta1)*MS1(4);
ydot(l+6:2*l+1) = k2*MS1(4:l-1)+a*M*S1(5:l)-(d+k2+delta2)*MS1(5:l);
ydot(2*l+2) = k2*MS1(l)-delta2*MS1(l+1);
ydot(2*l+3) = a*D*S1(2)-(d+k2+delta1)*DS1(1);
ydot(2*l+4) = k2*DS1(3)+a*D*S1(3)-(d+k2+delta1)*DS1(2);
ydot(2*l+5) = k2*DS1(4)+a*D*S1(4)-(d+k2+delta1)*DS1(3);
ydot(2*l+6:3*l+1) = k2*DS1(5:l)+a*D*S1(5:l)-(d+k2+delta2)*DS1(4:l-1);
ydot(3*l+2) = a*D*S1(l+1)-(d+k2+delta2)*DS1(l);
ydot(3*l+3) = Q2+d*MS2(1)-(a*M+delta1)*S2(1)+k2*DS2(1);
ydot(3*l+4) = d*(MS2(2)+DS2(1))-(a*M+a*D+delta1)*S2(2);
ydot(3*l+5) = d*(MS2(3)+DS2(2))-(a*M+a*D+delta1)*S2(3);
ydot(3*l+6) = d*(MS2(4)+DS2(3))-(a*M+a*D+delta1)*S2(4);
ydot(3*l+7:4*l+2) = d*(MS2(5:l)+DS2(4:l-1))-(a*M+a*D+delta2)*S2(5:l);
ydot(4*l+3) = d*DS2(l)-(a*D+delta2)*S2(l+1);
ydot(4*l+4) = a*M*S2(1)-(d+k1+delta1)*MS2(1);
ydot(4*l+5) = k1*MS2(1)+a*M*S2(2)-(d+k2+delta1)*MS2(2);
ydot(4*l+6) = k2*MS2(2)+a*M*S2(3)-(d+k2+delta1)*MS2(3);
ydot(4*l+7) = k2*MS2(3)+a*M*S2(4)-(d+k2+delta1)*MS2(4);
ydot(4*l+8:5*l+3) = k2*MS2(4:l-1)+a*M*S2(5:l)-(d+k2+delta2)*MS2(5:l);
ydot(5*l+4) = k2*MS2(l)-delta2*MS2(l+1);
ydot(5*l+5) = a*D*S2(2)-(d+k2+delta1)*DS2(1);
ydot(5*l+6) = k2*DS2(3)+a*D*S2(3)-(d+k2+delta1)*DS2(2);
ydot(5*l+7) = k2*DS2(4)+a*D*S2(4)-(d+k2+delta1)*DS2(3);
ydot(5*l+8:6*l+3) = k2*DS2(5:l)+a*D*S2(5:l)-(d+k2+delta2)*DS2(4:l-1);
ydot(6*l+4) = a*D*S2(l+1)-(d+k2+delta2)*DS2(l);
ydot(6*l+5) = (d+delta1)*(MS1(1)+MS1(2)+MS1(3)+MS1(4)+MS2(1)+MS2(2)+MS2(3)+...
    MS2(4))+(d+delta2)*(runsumMS1+runsumMS2)-...
    a*M*(firstrunsumS1+firstrunsumS2)+delta2*(MS1(l+1)+MS2(l+1));
ydot(6*l+6) = (d+delta1)*(DS1(1)+DS1(2)+DS1(3)+DS2(1)+DS2(2)+DS2(3))+...
    (d+delta2)*(runsumDS1+runsumDS2)-a*D*(secondrunsumS1+secondrunsumS2)+...
    k2*(DS1(1)+DS1(2)+DS2(1)+DS2(2));
end