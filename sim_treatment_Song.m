function [N,L,T,u,v,t] = sim_treatment_Song(D,TD,tf,Ts,N0,L0,T0,model)

TI=TD/10;

t=0:Ts:tf;
v=get_v(t,D,TD,TI);

[N,L,T,u,t] = sim_Song(tf,Ts,N0,L0,T0,0,v,model);

end