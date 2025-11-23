n = 9;
Q_size = (n+1)^2;
f = @(x,y) (n+1)*x+y+1;
len_d = 3;
len_b1 = 3;
len_b2 = 3;
len_a = n;
len_h = n;
len_l = n;
len_t = 8;
arr_d = linspace(0.01,0.99,len_d);
arr_b1 = linspace(0.01,0.99,len_b1);
arr_b2 = linspace(0.01,0.99,len_b2);
mfpt = zeros(len_d,len_b1,len_b2,len_a,len_h,len_l,len_t);
count = 0;
tt = zeros(18144,1);

for ind_d = 1:len_d
    d = arr_d(ind_d);
    
    for ind_b1 = 1:len_b1
        b1 = d + arr_b1(ind_b1);
        
        for ind_b2 = 1:len_b2
            b2 = d + arr_b2(ind_b2);
            
            for ind_a = 2:n-1
                a = ind_a;
                
                for ind_h = a+1:n
                    h = ind_h;
                    
                    for ind_l = 1:a-1
                    	l = ind_l;
                        arr_start = [f(h,h), f(h,l), f(l,l), f(l,h), f(h,h), f(l,h), f(l,l), f(h,l)];
                        arr_end = [f(h,l), f(l,l), f(l,h) f(h,h), f(l,h), f(l,l), f(h,l), f(h,h)];
                    
                        for ind_t = 1:len_t
                            Q = generator(b1,b2,d,n,Q_size);
                            Q(1,:) = [];
                            Q(:,1) = [];
                            count = count+1;
                            tt(count) = condest(Q);
                            p0 = zeros(Q_size,1);
                            p0(arr_start(ind_t)) = 1;
                            p0(1) = [];
                            tmp = -Q\p0;
                            mfpt(ind_d,ind_b1,ind_b2,ind_a,ind_h,ind_l,ind_t) = sum(tmp)-tmp(arr_end(ind_t)-1);
                        end
                        
                    end
                    
                end
                
            end
            
        end
               
    end
                
end

prob1 = mfpt(:,:,:,:,:,:,2)./(mfpt(:,:,:,:,:,:,2)+mfpt(:,:,:,:,:,:,8)); % A = Pr(HL --> HH)
prob2 = mfpt(:,:,:,:,:,:,6)./(mfpt(:,:,:,:,:,:,4)+mfpt(:,:,:,:,:,:,6)); % B = Pr(LH --> HH)
prob = prob1 + prob2 - prob1 .* prob2; % A union B = Pr(HL or LH --> HH) = A + B - A intersect B
save main.mat mfpt prob1 prob2 prob;

function Q = generator(b1,b2,d,n,Q_size)
f = @(x,y) (n+1)*x+y+1;
Qnz = 0;
tau = 6/n^2;
alph = 4/n;

for x = 1:1:n-1
    for y = 0:1:n
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x+1,y);
        Q_val(Qnz)=b1*x;
    end
end

for x = 2:1:n-1
    for y = 0:1:n
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x+1,y);
        Q_val(Qnz)=(alph/2)*x*(x-1);
    end
end

for x = 0:1:n
    for y = 1:1:n-1
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x,y+1);
        Q_val(Qnz)=b2*y;
    end
end

for x = 0:1:n
    for y = 2:1:n-1
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x,y+1);
        Q_val(Qnz)=(alph/2)*y*(y-1);
    end
end

for x = 0:1:n-1
    for y = 1:1:n
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x+1,y-1);
        Q_val(Qnz)=d*y;
    end
end

for x = 1:1:n
    for y = 0:1:n-1
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x-1,y+1);
        Q_val(Qnz)=d*x;
    end
end

for x = 3:1:n
    for y = 0:1:n
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x-1,y);
        Q_val(Qnz)=(tau/6)*x*(x-1)*(x-2);
    end
end

for x = 1:1:n
    for y = 0:1:n
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x-1,y);
        Q_val(Qnz)=x;
    end
end

for x = 0:1:n
    for y = 3:1:n
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x,y-1);
        Q_val(Qnz)=(tau/6)*y*(y-1)*(y-2);
    end
end

for x = 0:1:n
    for y = 1:1:n
        Qnz = Qnz + 1;
        Q_row(Qnz)=f(x,y);
        Q_col(Qnz)=f(x,y-1);
        Q_val(Qnz)=y;
    end
end

Q = sparse(Q_row,Q_col,Q_val,Q_size,Q_size,Qnz);
Q = Q';
Q = Q-diag(sum(Q));
end