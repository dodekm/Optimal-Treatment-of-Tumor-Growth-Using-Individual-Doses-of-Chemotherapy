function [c,c_eq,dc,dc_eq] = Lyap_con(x,x0,P)

c_eq=[];
dc_eq=[];

c=0.5*(x-x0)'*P*(x-x0)-1;
dc=P*(x-x0);

end