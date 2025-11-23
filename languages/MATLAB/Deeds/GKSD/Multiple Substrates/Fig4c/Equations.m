function ydot = Equations(t,y)
global Q1 Q2 delta1 delta2 Moff Mon Mcat Doff Don Dcat
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