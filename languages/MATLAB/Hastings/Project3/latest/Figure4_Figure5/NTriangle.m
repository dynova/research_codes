n = 6;
Q_size = (n+1)^3;
f = @(x,y,z) (n+1)^2*x+(n+1)*y+z+1;
len_d = 3;
len_b = 3;
len_t = 24;
tmp = 2:2:24;
tmp = [tmp;tmp-1];
tmp = tmp(:)';
arr_d = linspace(0.01,0.99,len_d);
arr_b = linspace(0.01,0.99,len_b);
mfpt = zeros(len_d,len_b,len_t);
a = 3;
h = n;
l = 1;

for ind_d = 1:len_d
    d = arr_d(ind_d);
    
    for ind_b = 1:len_b
        b = 2*d + arr_b(ind_b);
        
        arr_start = [f(h,h,h),f(h,h,l),f(h,h,l),f(h,l,l),f(h,l,l),f(l,l,l),f(l,l,l),f(l,l,h),f(l,l,h),f(l,h,h),f(l,h,h),f(l,h,l),f(h,h,h),f(h,l,h),f(h,h,h),f(l,h,h),f(h,h,l),f(l,h,l),f(h,l,l),f(h,l,h),f(l,l,l),f(l,h,l),f(l,l,h),f(h,l,h)];
        arr_end = arr_start(tmp);
        
        for ind_t = 1:len_t
            t = ind_t;
            initial = arr_start(t);
            final = arr_end(t);
            Q = generator(b,d,n,Q_size);
            Q(1,:) = [];
            Q(:,1) = [];
            p0 = zeros(Q_size,1);
            p0(initial) = 1;
            p0(1) = [];
            v = -Q\p0;
            mfpt(ind_d,ind_b,ind_t) = sum(v)-v(final-1);
            
        end
        
    end
    
end

r = 1./mfpt;
e12 = r(:,:,1);
e21 = r(:,:,2);
e23 = r(:,:,3);
e32 = r(:,:,4);
e34 = r(:,:,5);
e43 = r(:,:,6);
e45 = r(:,:,7);
e54 = r(:,:,8);
e56 = r(:,:,9);
e65 = r(:,:,10);
e67 = r(:,:,11);
e76 = r(:,:,12);
e18 = r(:,:,13);
e81 = r(:,:,14);
e16 = r(:,:,15);
e61 = r(:,:,16);
e27 = r(:,:,17);
e72 = r(:,:,18);
e38 = r(:,:,19);
e83 = r(:,:,20);
e47 = r(:,:,21);
e74 = r(:,:,22);
e58 = r(:,:,23);
e85 = r(:,:,24);
save NTriangle.mat e12 e21 e23 e32 e34 e43 e45 e54 e56 e65 e67 e76 e18 e81 e16 e61 e27 e72 e38 e83 e47 e74 e58 e85;

function Q = generator(b,d,n,Q_size)
f = @(x,y,z) (n+1)^2*x+(n+1)*y+z+1;
Qnz = 0;
tau = 6/n^2;
alph = 4/n;

for x = 1:1:n-1
    for y = 0:1:n
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x+1,y,z);
            Q_val(Qnz)=b*x;
        end
    end
end

for x = 2:1:n-1
    for y = 0:1:n
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x+1,y,z);
            Q_val(Qnz)=(alph/2)*x*(x-1);
        end
    end
end

for x = 0:1:n
    for y = 1:1:n-1
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y+1,z);
            Q_val(Qnz)=b*y;
        end
    end
end

for x = 0:1:n
    for y = 2:1:n-1
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y+1,z);
            Q_val(Qnz)=(alph/2)*y*(y-1);
        end
    end
end

for x = 0:1:n
    for y = 0:1:n
        for z = 1:1:n-1
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y,z+1);
            Q_val(Qnz)=b*z;
        end
    end
end

for x = 0:1:n
    for y = 0:1:n
        for z = 2:1:n-1
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y,z+1);
            Q_val(Qnz)=(alph/2)*z*(z-1);
        end
    end
end

for x = 3:1:n
    for y = 0:1:n
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x-1,y,z);
            Q_val(Qnz)=(tau/6)*x*(x-1)*(x-2);
        end
    end
end

for x = 0:1:n
    for y = 3:1:n
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y-1,z);
            Q_val(Qnz)=(tau/6)*y*(y-1)*(y-2);
        end
    end
end

for x = 0:1:n
    for y = 0:1:n
        for z = 3:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y,z-1);
            Q_val(Qnz)=(tau/6)*z*(z-1)*(z-2);
        end
    end
end

for x = 1:1:n
    for y = 0:1:n
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x-1,y,z);
            Q_val(Qnz)=x;
        end
    end
end

for x = 0:1:n
    for y = 1:1:n
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y-1,z);
            Q_val(Qnz)=y;
        end
    end
end

for x = 0:1:n
    for y = 0:1:n
        for z = 1:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y,z-1);
            Q_val(Qnz)=z;
        end
    end
end

for x = 0:1:n-1
    for y = 1:1:n
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x+1,y-1,z);
            Q_val(Qnz)=d*y;
        end
    end
end

for x = 1:1:n
    for y = 0:1:n-1
        for z = 0:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x-1,y+1,z);
            Q_val(Qnz)=d*x;
        end
    end
end

for x = 0:1:n
    for y = 1:1:n
        for z = 0:1:n-1
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y-1,z+1);
            Q_val(Qnz)=d*y;
        end
    end
end

for x = 0:1:n
    for y = 0:1:n-1
        for z = 1:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x,y+1,z-1);
            Q_val(Qnz)=d*z;
        end
    end
end

for x = 1:1:n
    for y = 0:1:n
        for z = 0:1:n-1
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x-1,y,z+1);
            Q_val(Qnz)=d*x;
        end
    end
end

for x = 0:1:n-1
    for y = 0:1:n
        for z = 1:1:n
            Qnz = Qnz + 1;
            Q_row(Qnz)=f(x,y,z);
            Q_col(Qnz)=f(x+1,y,z-1);
            Q_val(Qnz)=d*z;
        end
    end
end

Q = sparse(Q_row,Q_col,Q_val,Q_size,Q_size,Qnz);
Q = Q';
Q = Q-diag(sum(Q));
end