function [Nr,Lr,Tr,u0] = get_ss_Song(model,v0)

[a,b,c,d,r,mu_,alpha_1,alpha_2,beta_1,beta_2,omega,k_N,k_L,k_T] = parameters_Song(model);

u0=v0/omega;

syms Nr Lr Tr

eqs=[Nr*(a*(1-b*Nr)-alpha_1*Tr-k_N*u0)==0;
    r*Nr*Tr-(mu_+beta_1*Tr+k_L*u0)*Lr==0;
    Tr*(c*(1-d*Tr)-alpha_2*Nr-beta_2*Lr-k_T*u0)==0];

S=solve(eqs,[Nr,Lr,Tr]);
Nr=double(vpa(S.Nr));
Lr=double(vpa(S.Lr));
Tr=double(vpa(S.Tr));

idx=(Nr>=0 & Lr>=0 & Tr>=0 & imag(Nr)==0 & imag(Lr)==0 & imag(Tr)==0);
Nr=Nr(idx);
Lr=Lr(idx);
Tr=Tr(idx);

end