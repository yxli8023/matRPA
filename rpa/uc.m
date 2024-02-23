function ucc = uc(nq)
global chi;
I=eye(16);
U=zeros(16);
for u=0.0:0.1:10.0
    U(1,1)=u;U(6,6)=u;U(11,11)=u;U(16,16)=u;
    for iq=1:nq
        chi0=chi{iq};
        chis=I-U*chi0;
        [~,E1]=eig(chis);
        e=[min(real(diag(E1)))];   
        if min(e)<0
            break
        end
    end
    if min(e)<0
        break
    end
end 
ucc=u;
end