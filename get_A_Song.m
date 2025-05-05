function A = get_A_Song(x,model)

[a,b,c,d,r,mu_,alpha_1,alpha_2,beta_1,beta_2,omega,k_N,k_L,k_T] = parameters_Song(model);

N=x(1);
L=x(2);
T=x(3);
u=x(4);

A=[a*(1-2*b*N)-alpha_1*T-k_N*u, 0 , -alpha_1*N, -k_N*N; 
    r*T , -(mu_+beta_1*T+k_L*u) , r*N-beta_1*L, -k_L*L ;
    -alpha_2*T , -beta_2*T , c*(1-2*d*T)-alpha_2*N-beta_2*L-k_T*u, -k_T*T ;
    0,0,0, -omega ];

end