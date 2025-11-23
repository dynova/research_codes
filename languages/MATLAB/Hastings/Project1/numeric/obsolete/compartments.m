%% Preliminaries
N = 25;
Q_size = (N+1)^2;
idx = @(x,y) (N+1)*x+y+1;
[array_b1, array_b2, array_d] = deal(linspace(0.01,0.99,2)*N);
array_t = logspace(-2,6,1000);
len_b1 = length(array_b1);
len_b2 = length(array_b2);
len_d = length(array_d);
len_t = length(array_t);
len_HL = 8;
len_init = 4;
array_init = eye(len_init,len_HL/2);
L1 = 1;
L2 = 1;
H1 = floor(1+sqrt(0.5*N));
H2 = floor(1+sqrt(0.5*N));
mfpt = zeros(len_d,len_b1,len_b2,len_HL);
state_prob = zeros(len_d,len_b1,len_b2,len_t,len_init,len_HL/2);
r = zeros(len_d,len_b1,len_b2,len_HL);
cmp = cell(len_d,len_b1,len_b2);

%% Begin loops
for ind_d = 1:len_d
    d = array_d(ind_d);
    
    for ind_b1 = 1:len_b1
        b1 = array_b1(ind_b1);
        
        for ind_b2 = 1:len_b2
            b2 = array_b2(ind_b2);
            
            array_start = ...
                {idx(H1,H2), idx(H1,H2), idx(H1,L2), idx(L1,H2),...
                idx(L1,H2) idx(H1,L2) idx(L1,L2) idx(L1,L2)
                };
            
            array_end = ...
                {idx(L1,H2), idx(H1,L2), idx(L1,L2), idx(L1,L2),...
                idx(H1,H2) idx(H1,H2) idx(H1,L2) idx(L1,H2)
                };
            
            Q = gen(b1,b2,d,N,Q_size);
            
            for ind_HL = 1:len_HL
                S = full(Q);
                S([1,array_end{ind_HL}],:) = [];
                S(:,[1,array_end{ind_HL}]) = [];
                S = sparse(S);
                p_0 = zeros(Q_size,1);
                p_0(array_start{ind_HL}) = 1;
                p_0([1,array_end{ind_HL}]) = [];
                p_0 = sparse(p_0);
                mfpt(ind_d,ind_b1,ind_b2,ind_HL) = ...
                    sum(-S\p_0);
                r(ind_d,ind_b1,ind_b2,ind_HL) = 1/mfpt(ind_d,ind_b1,ind_b2,ind_HL);
                cmp{ind_d,ind_b1,ind_b2} = ...
                    [-r(1)-r(2) r(6)       r(5)        0;         ...
                    r(2)      -r(3)-r(6) 0            r(7);       ...
                    r(1)       0          -r(4)-r(5)  r(8);       ...
                    0           r(3)       r(4)       -r(7)-r(8); ...
                    ];
                
                for ind_init = 1:len_init
                    init = array_init(ind_init,:);
                    
                    for ind_t = 1:len_t
                        t = array_t(ind_t);
                        t4 = expm(cmp{ind_d,ind_b1,ind_b2}*t)*init';
                        state_prob(ind_d,ind_b1,ind_b2,ind_t,ind_init,:) = t4./sum(t4);
                    end
                end
            end
        end
    end
end

save compartments.mat;
save mma.mat cmp r;

%% Subroutine 'gen'
function Q = gen(b1,b2,d,N,Q_size)
idx = @(x,y) (N+1)*x+y+1;
Qnz = 0;
for x = 1:1:N-1
    for y = 0:1:N
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x+1,y);
        Q_val(Qnz)=b1*x;
    end
end

for x = 2:1:N-1
    for y = 0:1:N
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x+1,y);
        Q_val(Qnz)=2*x*(x-1);
    end
end

for x = 0:1:N
    for y = 1:1:N-1
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x,y+1);
        Q_val(Qnz)=b2*y;
    end
end

for x = 0:1:N
    for y = 2:1:N-1
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x,y+1);
        Q_val(Qnz)=2*y*(y-1);
    end
end

for x = 0:1:N-1
    for y = 1:1:N
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x+1,y-1);
        Q_val(Qnz)=d*y;
    end
end

for x = 1:1:N
    for y = 0:1:N-1
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x-1,y+1);
        Q_val(Qnz)=d*x;
    end
end

for x = 3:1:N
    for y = 0:1:N
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x-1,y);
        Q_val(Qnz)=x*(x-1)*(x-2);
    end
end

for x = 1:1:N
    for y = 0:1:N
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x-1,y);
        Q_val(Qnz)=x;
    end
end

for x = 0:1:N
    for y = 3:1:N
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x,y-1);
        Q_val(Qnz)=y*(y-1)*(y-2);
    end
end

for x = 0:1:N
    for y = 1:1:N
        Qnz = Qnz + 1;
        Q_row(Qnz)=idx(x,y);
        Q_col(Qnz)=idx(x,y-1);
        Q_val(Qnz)=y;
    end
end

Q=sparse(Q_row,Q_col,Q_val,Q_size,Q_size,Qnz);
Q = Q';
for i = 1:Q_size
    Q(i,i) = Q(i,i)-sum(Q(:,i));
end
end