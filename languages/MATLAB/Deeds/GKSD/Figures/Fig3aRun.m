n = 20;
rng = 1e6;
Q = 2e-2;
Qfinal = Q*rng;
inc = 10^(log10(rng)/(n-1));
len = 100;
rvec = logspace(0,6,len);
ind = 1;
r1 = 2e-5;
y = zeros(len,n);
while ind <= n
    S00 = Q/r1;
    y(:,ind) = Fig3aScript(Q,S00);
    Q = Q*inc;
    ind = ind + 1;
end
save 3a