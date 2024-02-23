function A = matr(tpara,kx,ky,kz)
t1 = tpara(1); t2 = tpara(2); t3 = tpara(3); t4 = tpara(4);
s=zeros(4,4);
s(1,2) = -1*t1;
s(1,3) = -1*t2*exp((1i)*ky)-t3;
s(1,4) = -1*t1;
s(2,3) = -1*t1;
s(2,4) = -1*t2*exp((1i)*kx)-t3;
s(3,4) = -1*t1;
s(1,1) = -1*t4*cos(kz);
s(2,2) = -1*t4*cos(kz);
s(3,3) = -1*t4*cos(kz);
s(4,4) = -1*t4*cos(kz);
A = s+s';
end

function [A,B]= eigs(tpara,nk,nkz)
eng=cell(4,1);
mat=cell(16,1);
for ix=1:nk
    kx=2*pi*(ix-1)/(nk);
    for iy=1:nk
        ky=2*pi*(iy-1)/(nk);
        for iz=1:nkz
            kz=2*pi*(iz-1)/(nkz);
            [M,E]=eig(matr(tpara,kx,ky,kz));
            for j=1:4
                eng{j,1}(ix,iy,iz) = E(j,j);
            end
            for j=1:16
                mat{j,1}(ix,iy,iz)=M(j);
            end
        end
    end
end
A=eng;
B=mat;
end


