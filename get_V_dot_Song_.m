function [V_dot,dV_dot] = get_V_dot_Song_(x,model,x0,P)

[a,b,c,d,r,mu_,alpha_1,alpha_2,beta_1,beta_2,omega,k_N,k_L,k_T] = parameters_Song(model);

dot_x=get_dX_Song(model,[x;0],0);
dot_x=dot_x(1:3);

x=x(1:3);

%V_dot=0.5*((x-x0)'*P*dot_x)+(dot_x'*P*(x-x0));

V_dot=((x-x0)'*P*dot_x);

A=get_A_Song([x;0],model);
A=A(1:3,1:3);
dV_dot=P*dot_x+A'*P*(x-x0);


end