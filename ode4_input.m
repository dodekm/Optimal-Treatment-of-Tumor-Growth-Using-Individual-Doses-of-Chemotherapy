function X = ode4_input(odefun,t,X0,U)
 
h = diff(t);

n = length(X0);
N = length(t);
X = zeros(n,N);
F = zeros(n,4);

X(:,1) = X0;
for i = 2:N
  hi = h(i-1);
  xi = X(:,i-1);
  F(:,1) = feval(odefun,xi,U(:,i));
  F(:,2) = feval(odefun,xi+0.5*hi*F(:,1),U(:,i));
  F(:,3) = feval(odefun,xi+0.5*hi*F(:,2),U(:,i));  
  F(:,4) = feval(odefun,xi+hi*F(:,3),U(:,i));
  X(:,i) = xi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
end

