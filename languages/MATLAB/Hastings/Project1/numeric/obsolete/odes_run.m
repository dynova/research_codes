global b1 b2 d
k = 50;
[N1, N2] = deal(k);
b1 = 0.99*N1;
b2 = 0.99*N2; 
d = 1;
y(1,:) = run(N1,N2);
      
function y = run(N1,N2)
[~,x] = ode15s(@odes,linspace(0,1e0,1e4),[N1 N2]);
y = x(end,:);
end

function ydot = odes(~,y)
global b1 b2 d
X1 = y(1);
X2 = y(2);
ydot = zeros(2,1);
ydot(1) = (b1-d-1).*X1 + 2*X1.^2 - X1.^3 + d.*X2;
ydot(2) = (b2-d-1).*X2 + 2*X2.^2 - X2.^3 + d.*X1;
end