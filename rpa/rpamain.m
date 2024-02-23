global chi;
global eng;
global mat;
global mu;
tpara=[2.686,2.985,0.574,0.259];
nk=10;
nkz=10;
fill=1.10;
u=4.5;
[eng,mat,mu]=eigs(tpara,nk,nkz,fill);
chi=chi2(nk,nkz);
uc=uc(nk*nk*nkz);
[chis,chic]=fermi(nk,nkz,u);

%归一化与对角化
[Lamdas,Deltas]=eig(chis);
[Lamdat,Deltac]=eig(chic);
save('data.mat','chi','chic','chis','Deltas','Deltac','uc','mu');
%Deltas=sort(real(diag(Deltas)),'descend');
%Deltac=sort(real(diag(Deltac)),'descend');