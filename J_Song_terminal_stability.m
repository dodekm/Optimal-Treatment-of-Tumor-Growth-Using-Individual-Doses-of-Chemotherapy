function [J] = J_Song_terminal_stability(D,TD,tf,Ts,w_N,w_L,w_T,w_u,w_s,P,Nr,Lr,Tr,N0,L0,T0,model)

[N,L,T,u] = sim_treatment_Song(D,TD,tf,Ts,N0,L0,T0,model);

J_N=sum((N-Nr).^2);
J_L=sum((L-Lr).^2);
J_T=sum((T-Tr).^2);
J_u=sum(u.^2);

J_s=max(Lyap_con([N(end);L(end);T(end)],[Nr;Lr;Tr],P),0);
J=J_N*w_N+J_L*w_L+J_T*w_T+J_u*w_u+J_s*w_s;

end