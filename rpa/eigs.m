function [A,B,mu]= eigs(tpara,nk,nkz,fill)
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
eng2=[eng{1},eng{2},eng{3},eng{4}];
eng2=sort(eng2(:));
mu=eng2(round(length(eng2)*fill/2));
A=eng;
B=mat;
end