function out = acf(x)
n = length(x); m = mean(x); num = 0; denom = 0;
for i = 1:n-1
    num = num+(x(i)-m)*(x(i+1)-m);
    denom = denom+(x(i)-m)*(x(i)-m);
end
denom = denom+(x(n)-m)*(x(n)-m);
out = num/denom;
end