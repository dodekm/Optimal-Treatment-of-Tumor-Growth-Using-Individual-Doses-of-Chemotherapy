function [N,L,T,u,t] = sim_Song(tf,Ts,N0,L0,T0,u0,v,model)

t=0:Ts:tf;

dx_fcn=@(x,u) get_dX_Song(model,x,u);
x0=[N0;L0;T0;u0];

x = ode4_input(dx_fcn,t,x0,v)';

N=x(:,1);
L=x(:,2);
T=x(:,3);
u=x(:,4);

end