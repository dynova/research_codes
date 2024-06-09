function ydot = EquationsDDT2(t,y)
global Q1 Q2 r1 r2 d a k l
S1 = y(1:l+1);
MS1 = y(l+2:2*l+1);
DS1 = y(2*l+2:3*l+1);
S2 = y(3*l+2:4*l+2);
MS2 = y(4*l+3:5*l+2);
DS2 = y(5*l+3:6*l+2);
M = y(6*l+3);
D = y(6*l+4);

ydot = zeros(6*l+4,1);
firstrunsumDS1 = sum(DS1(1:l));
firstrunsumS1 = sum(S1(1:l));
secondrunsumS1 = sum(S1(2:l+1));
firstrunsumDS2 = sum(DS2(1:l));
firstrunsumS2 = sum(S2(1:l));
secondrunsumS2 = sum(S2(2:l+1));
secondrunsumDS1 = sum(DS1(4:l));
secondrunsumDS2 = sum(DS2(4:l));
runsumMS1 = sum(MS1(5:l));
runsumMS2 = sum(MS2(5:l));

ydot(1) = Q1-r1*S1(1)+d*MS1(1)-a*M*S1(1)+k*firstrunsumDS1;
ydot(2) = -S1(2)*(r1+a*M+a*D)+k*MS1(1)+d*(DS1(1)+MS1(2));
ydot(3) = -S1(3)*(r1+a*M+a*D)+k*MS1(2)+d*(DS1(2)+MS1(3));
ydot(4) = -S1(4)*(r1+a*M+a*D)+k*MS1(3)+d*(DS1(3)+MS1(4));
ydot(5:l) = -S1(5:l)*(r2+a*M+a*D)+k*MS1(4:l-1)+d*(DS1(4:l-1)+MS1(5:l));
ydot(l+1) = k*MS1(l)-(a*D+r2)*S1(l+1)+d*DS1(l);
ydot(l+2) = a*M*S1(1)-(d+k+r1)*MS1(1);
ydot(l+3) = a*M*S1(2)-(d+k+r1)*MS1(2);
ydot(l+4) = a*M*S1(3)-(d+k+r1)*MS1(3);
ydot(l+5) = a*M*S1(4)-(d+k+r1)*MS1(4);
ydot(l+6:2*l+1) = a*M*S1(5:l)-(d+k+r2)*MS1(5:l);
ydot(2*l+2) = a*D*S1(2)-(d+k+r1)*DS1(1);
ydot(2*l+3) = a*D*S1(3)-(d+k+r1)*DS1(2);
ydot(2*l+4) = a*D*S1(4)-(d+k+r1)*DS1(3);
ydot(2*l+5:3*l+1) = a*D*S1(5:l+1)-(d+k+r2)*DS1(4:l);
ydot(3*l+2) = Q2-r1*S2(1)+d*MS2(1)-a*M*S2(1)+k*firstrunsumDS2;
ydot(3*l+3) = -S2(2)*(r1+a*M+a*D)+k*MS2(1)+d*(DS2(1)+MS2(2));
ydot(3*l+4) = -S2(3)*(r1+a*M+a*D)+k*MS2(2)+d*(DS2(2)+MS2(3));
ydot(3*l+5) = -S2(4)*(r1+a*M+a*D)+k*MS2(3)+d*(DS2(3)+MS2(4));
ydot(3*l+6:4*l+1) = -S2(5:l)*(r2+a*M+a*D)+k*MS2(4:l-1)+d*(DS2(4:l-1)+MS2(5:l));
ydot(4*l+2) = k*MS2(l)-(a*D+r2)*S2(l+1) + d*DS2(l);
ydot(4*l+3) = a*M*S2(1)-(d+k+r1)*MS2(1);
ydot(4*l+4) = a*M*S2(2)-(d+k+r1)*MS2(2);
ydot(4*l+5) = a*M*S2(3)-(d+k+r1)*MS2(3);
ydot(4*l+6) = a*M*S2(4)-(d+k+r1)*MS2(4);
ydot(4*l+7:5*l+2) = a*M*S2(5:l)-(d+k+r2)*MS2(5:l);
ydot(5*l+3) = a*D*S2(2)-(d+k+r1)*DS2(1);
ydot(5*l+4) = a*D*S2(3)-(d+k+r1)*DS2(2);
ydot(5*l+5) = a*D*S2(4)-(d+k+r1)*DS2(3);
ydot(5*l+6:6*l+2) = a*D*S2(5:l+1)-(d+k+r2)*DS2(4:l);
ydot(6*l+3) = (d+k+r1)*(MS1(1)+MS1(2)+MS1(3)+MS1(4)+MS2(1)+MS2(2)+MS2(3)+MS2(4))+ ...
    + (d+k+r2)*(runsumMS1+runsumMS2)-a*M*(firstrunsumS1+firstrunsumS2);
ydot(6*l+4) = (d+k+r1)*(DS1(1)+DS1(2)+DS1(3)+DS2(1)+DS2(2)+DS2(3)) + ...
    + (d+k+r2)*(secondrunsumDS1+secondrunsumDS2)-a*D*(secondrunsumS1+secondrunsumS2);
end