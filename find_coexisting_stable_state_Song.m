function [Nr,Lr,Tr] = find_coexisting_stable_state_Song(model,Nr,Lr,Tr,ur)

i=(Nr>=1e3 & Lr>=1e3 & Tr>=1e3 & imag(Nr)==0 & imag(Lr)==0 & imag(Tr)==0);
Nr=Nr(i);
Lr=Lr(i);
Tr=Tr(i);

for i=1:length(Nr)
    
    A = get_A_Song([Nr(i);Lr(i);Tr(i);ur],model);

    if sum(real(eig(A))<0)==4
        break;
    end    

end


Nr=Nr(i);
Lr=Lr(i);
Tr=Tr(i);

end
