function dx = get_dX_Song(model,x,v)

[a,b,c,d,r,mu_,alpha_1,alpha_2,beta_1,beta_2,omega,k_N,k_L,k_T] = parameters_Song(model);

N=x(1);
L=x(2);
T=x(3);
u=x(4);

dN=N*(a*(1-b*N)-alpha_1*T-k_N*u);
dL=r*N*T-L*(mu_+beta_1*T+k_L*u);
dT=T*(c*(1-d*T)-alpha_2*N-beta_2*L-k_T*u);

% dN=a*N*(1-b*N)-alpha_1*N*T-k_N*u*N;
% dL=r*N*T-mu_*L-beta_1*L*T-k_L*u*L;
% dT=c*T*(1-d*T)-alpha_2*N*T-beta_2*L*T-k_T*u*T;
du=v-omega*u;

dx=[dN;dL;dT;du];

end