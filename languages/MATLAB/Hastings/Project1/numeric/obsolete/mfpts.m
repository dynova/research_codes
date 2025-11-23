N = 10;
Q_size = (N+1)^2;
idx = @(x,y) (N+1)*x+y+1;
[array_b1, array_b2] = deal(linspace(0.01,0.99,2)*N);
array_d = linspace(0.01,0.99,2)*N;
len_b1 = length(array_b1);
len_b2 = length(array_b2);
len_d = length(array_d);
array_Q = cell(len_d,len_b1,len_b2);
len_t = 1;
L1 = 1;
L2 = 1;
H1 = floor(1+sqrt(0.5*N));
H2 = floor(1+sqrt(0.5*N));
array_start = [idx(H1,H2), idx(H1,H2), idx(H1,L2), idx(L1,H2)];
array_end = [idx(L1,H2), idx(H1,L2), idx(L1,L2), idx(L1,L2)];
mfpt = zeros(len_d,len_b1,len_b2,len_t);
array_S = cell(len_d,len_b1,len_b2,len_t);
array_p0 = cell(len_d,len_b1,len_b2,len_t);

%% Begin loops
for ind_d = 1:len_d
    d = array_d(ind_d);
    
    for ind_b1 = 1:len_b1
        b1 = array_b1(ind_b1);
        
        for ind_b2 = 1:len_b2
            b2 = array_b2(ind_b2);
            Q = gen(b1,b2,d,N,Q_size);
                        
            for ind_t = 1:len_t
                S = full(Q);
                S([1,array_end(ind_t)],:) = [];
                S(:,[1,array_end(ind_t)]) = [];
                S = sparse(S);
                p0 = zeros(Q_size,1);
                p0(array_start(ind_t)) = 1;
                p0([1,array_end(ind_t)]) = [];
                p0 = sparse(p0);
                mfpt(ind_d,ind_b1,ind_b2,ind_t) = sum(-S\p0);
                array_S{ind_d,ind_b1,ind_b2,ind_t} = S;
                array_p0{ind_d,ind_b1,ind_b2,ind_t} = p0;
            end
        end
    end
end

% save mfpts.mat mfpt array_S array_p0;

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