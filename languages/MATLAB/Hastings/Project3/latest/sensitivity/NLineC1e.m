function val = NLineC1e(n)
len_g = 5;
g_vec = linspace(0.75,0.95,len_g);
pp = 10;
T = Inf;
S = simplexgrid(8,pp,1);
load NLine.mat e12 e21 e23 e32 e34 e43 e45 e54 e56 e65 e67 e76 e18 e81 e16 e61 e27 e72 e38 e83 e47 e74 e58 e85;
r1 = e12;
r2 = e21;
r3 = e23;
r4 = e32;
r5 = e34;
r6 = e43;
r7 = e45;
r8 = e54;
r9 = e56;
r10 = e65;
r11 = e67;
r12 = e76;
r13 = e18;
r14 = e81;
r15 = e16;
r16 = e61;
r17 = e27;
r18 = e72;
r19 = e38;
r20 = e83;
r21 = e47;
r22 = e74;
r23 = e58;
r24 = e85;
c1 = 2;
c2 = 25;
len_d = 3;
len_b = 3;
len_a = n;
len_h = n;
len_l = n;
val = zeros(len_g,len_d,len_b,len_a,len_h,len_l);
count = 0;
Xvals = [repelem(1:4,1,8)' repmat(1:8,1,4)'];
P1 = [0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 1, 0, 0, 0];
P3 = [0, 1, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 1, 0, 0, 0, 0, 0];
P4 = [0, 0, 0, 0, 0, 0, 0, 1;
    0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1];
R = [0 c1*c2 0 0 c2 c2+c1*c2 c2 c2 2*c2 2*c2+c1*c2 2*c2 2*c2 3*c2 3*c2+c1*c2 3*c2 3*c2 2*c2 2*c2+c1*c2 2*c2 2*c2 c2 c2+c1*c2 c2 c2 2*c2 2*c2+c1*c2 2*c2 2*c2 c2 c2+c1*c2 c2 c2]-c1*c2;
opt1 = struct('Qtype',0,'Rtype',2);
opt2 = struct('maxit',600,'nochangelim',500,'prtiters',0,'print',0);
p = 1;
q = 0;
Q1 = [p q;
    p q;
    p q;
    q p;
    q p;
    q p;
    q p;
    p q];
Q = repmat(Q1',1,4);

for ind_d = 1:len_d

    for ind_b = 1:len_b

        for ind_a = 2:n-1
            a = ind_a;

            for ind_h = a+1:n
                h = ind_h;

                for ind_l = 1:a-1
                    l = ind_l;

                    e12 = r1(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e21 = r2(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e23 = r3(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e32 = r4(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e34 = r5(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e43 = r6(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e45 = r7(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e54 = r8(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e56 = r9(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e65 = r10(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e67 = r11(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e76 = r12(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e18 = r13(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e81 = r14(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e16 = r15(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e61 = r16(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e27 = r17(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e72 = r18(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e38 = r19(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e83 = r20(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e47 = r21(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e74 = r22(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e58 = r23(ind_d,ind_b,ind_a,ind_h,ind_l);
                    e85 = r24(ind_d,ind_b,ind_a,ind_h,ind_l);
                    P2 = [0, e12/(e12+e16+e18), 0, 0, 0, e16/(e12+e16+e18), 0, e18/(e12+e16+e18);
                        e21/(e21+e23+e27), 0, e23/(e21+e23+e27), 0, 0, 0, e27/(e21+e23+e27), 0;
                        0, e32/(e32+e34+e38), 0, e34/(e32+e34+e38), 0, 0, 0, e38/(e32+e34+e38);
                        0, 0, e43/(e43+e45+e47), 0, e45/(e43+e45+e47), 0, e47/(e43+e45+e47), 0;
                        0, 0, 0, e54/(e54+e56+e58), 0, e56/(e54+e56+e58), 0, e58/(e54+e56+e58);
                        e61/(e61+e65+e67), 0, 0, 0, e65/(e61+e65+e67), 0, e67/(e61+e65+e67), 0;
                        0, e72/(e72+e74+e76), 0, e74/(e72+e74+e76), 0, e76/(e72+e74+e76), 0, 0;
                        e81/(e81+e83+e85), 0, e83/(e81+e83+e85), 0, e85/(e81+e83+e85), 0, 0, 0];
                    P = [P1' P2' P3' P4'];
                    [b,Pb,Rb] = pomdp(pp,P,Q,R,opt1);

                    for ind_g = 1:len_g
                        gma = g_vec(ind_g);
                        model = struct('P',Pb,'R',Rb,'discount',gma,'T',T);
                        results = mdpsolve(model,opt2);
                        f = results.v;
                        val(ind_g,ind_d,ind_b,ind_a,ind_h,ind_l) = f(S(:,4)==1);
                        count = count + 1;
%                         disp(count)
                    end

                end

            end

        end

    end

end

% save NLineC1e.mat val;