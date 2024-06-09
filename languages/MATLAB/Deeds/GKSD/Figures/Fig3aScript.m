function Stot = Fig3aScript(Q,S00)
global a d k r1 r2

E3vec = logspace(-1,5,100);
lenE3 = length(E3vec);
S0plot = zeros(lenE3,1);
S1plot = zeros(lenE3,1);
E3S0plot = zeros(lenE3,1);
BS1plot = zeros(lenE3,1);
r1 = 2e-5;
r2 = 2e-4;
a = 1e-4;
d = 1e-3;
k = 0.999;
S10 = 0;
E3S00 = 0;
BS10 = 0;
B0 = 0.1;
tf = 1e10;
time = 0:tf/1E3:tf;

for E3loop = 1:lenE3    
    E30 = E3vec(E3loop);
    y0 = [S00 S10 E3S00 BS10 E30 B0];
    [t,y] = ode23s(@GKSynDegEqns,time,y0);
    S0plot(E3loop,1) = y(lenE3,1);
    S1plot(E3loop,1) = y(lenE3,2);
    E3S0plot(E3loop,1) = y(lenE3,3);
    BS1plot(E3loop,1) = y(lenE3,4);
    Stot = S0plot + S1plot + E3S0plot + BS1plot;
end
    function ydot = GKSynDegEqns(t,y)
        S0 = y(1);
        S1 = y(2);
        E3S0 = y(3);
        BS1 = y(4);
        E3 = y(5);
        B = y(6);
        ydot = zeros(6,1);
        ydot(1) = Q-(r1+a*E3)*S0 + d*E3S0 + k*BS1;
        ydot(2) = k*E3S0 + d*BS1 - (a*B+r2)*S1;
        ydot(3) = a*E3*S0 - (d+k+r1)*E3S0;
        ydot(4) = a*B*S1 - (d+k+r2)*BS1;
        ydot(5) = -a*E3*S0 + (d+k+r1)*E3S0;
        ydot(6) = -a*B*S1 + (d+k+r2)*BS1;
    end
end