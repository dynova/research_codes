load Fig4.mat;
ind = 3;
period = 10;
y = squeeze(x2(ind,:,1:period:end));
x = 0:time_max*period/num_steps:time_max;
n = length(x);
z = zeros(num_sim,n);
h = 50;
wfun = @(h,mu,x) (2*pi*h^2)^(-1/2)*exp(-0.5*((x-mu)/h).^2);

for ind = 1:n
    w = wfun(h,x(ind),x);
    for ind_sim = 1:num_sim
        z(ind_sim,ind) = sum(w.*y(ind_sim,:))/sum(w);
    end
end

res = y - z;
n1 = size(res,1);
n2 = size(res,2);
w = floor(n2/2);

temp1 = zeros(n1,n2);
for i = 1:n1
    for j = w:n2
        temp1(i,j) = acf(res(i,(j-w+1:j)));
    end
end
acf1 = temp1(:,w:n2);
temp2 = quantile(acf1,[0.025,0.5,0.975]);
acf1_lb = temp2(1,:);
acf1_mdn = temp2(2,:);
acf1_ub = temp2(3,:);
acf1_t = median(corr((w:n2)',acf1','type','Kendall'));

temp1 = zeros(n1,n2);
for i = 1:n1
    for j = w:n2
        temp1(i,j) = var(res(i,(j-w+1:j)));
    end
end
vr = temp1(:,w:n2);
temp2 = quantile(vr,[0.025,0.5,0.975]);
vr_lb = temp2(1,:);
vr_mdn = temp2(2,:);
vr_ub = temp2(3,:);
vr_t = median(corr((w:n2)',vr','type','Kendall'));

temp1 = zeros(n1,n2);
temp2 = zeros(n1,n2);
for i = 1:n1
    for j = w:n2
        temp1(i,j) = std(y(i,(j-w+1:j)));
        temp2(i,j) = mean(y(i,(j-w+1:j)));
    end
end
cv = temp1(:,w:n2)./temp2(:,w:n2);
temp2 = quantile(cv,[0.025,0.5,0.975]);
cv_lb = temp2(1,:);
cv_mdn = temp2(2,:);
cv_ub = temp2(3,:);
cv_t = median(corr((w:n2)',cv','type','Kendall'));

save Fig8_23.mat acf1_lb acf1_mdn acf1_ub acf1_t...
                 vr_lb vr_mdn vr_ub vr_t...
                 cv_lb cv_mdn cv_ub cv_t;