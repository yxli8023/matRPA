function [chis,chic] = fermi(nk,nkz,u)
global chi;
global eng;
global mat;
global mu;
nq = nk * nk * nkz;
vq2 = [];
vq3 = [];
eng2 = zeros(nk,nk,nkz);
eng3 = zeros(nk,nk,nkz);
for iq = 1:nq
    iqx = mod(iq-1,nk) + 1;
    iqy = floor( (mod(iq-1,nk * nk))/nk ) + 1;
    iqz = floor( (iq-1)/(nk * nk) ) + 1;
    eng2(iqx,iqy,iqz) = eng{2}(iq);
    eng3(iqx,iqy,iqz) = eng{3}(iq);
    dqy = min(nk,iqy + 1);
    dqx = min(nk,iqx + 1);
    nq1 = (iqz - 1) * nk * nk + (dqy-1) * nk + iqx;
    nq2 = (iqz - 1) * nk * nk + (iqy-1) * nk+dqx;
    if (eng{2}(iq)-mu) * (eng{2}(nq1)-mu)<0
        if abs(eng{2}(iq))<abs(eng{2}(nq1))
            vq2 = [vq2,iq];
        else
            vq2 = [vq2,nq1];
        end
    end
    if (eng{3}(iq)-mu) * (eng{3}(nq1)-mu)<0
        if abs(eng{3}(iq))<abs(eng{3}(nq1))
            vq3 = [vq3,iq];
        else
            vq3 = [vq3,nq1];
        end
    end
    if (eng{2}(iq)-mu) * (eng{2}(nq2)-mu)<0
        if abs(eng{2}(iq))<abs(eng{2}(nq2))
            vq2 = [vq2,iq];
        else
            vq2 = [vq2,nq2];
        end
    end
    if (eng{3}(iq)-mu) * (eng{3}(nq2)-mu)<0
        if abs(eng{3}(iq))<abs(eng{3}(nq2))
        vq3 = [vq3,iq];
        else
        vq3 = [vq3,nq2];
        end
    end
    vq2 = unique(vq2);
    vq3 = unique(vq3);
end
[Fx,Fy,Fz]  =  gradient(real(eng2));
eng2 = sqrt(Fx.^2 * nk * nk + Fy.^2 * nk * nk + Fz.^2 * nkz * nkz);
%eng2 = sqrt(Fx.^2 * nk * nk + Fy.^2 * nk * nk);
[Fx,Fy,Fz]  =  gradient(real(eng3));
eng3 = sqrt(Fx.^2 * nk * nk + Fy.^2 * nk * nk + Fz.^2 * nkz * nkz);
%eng3 = sqrt(Fx.^2 * nk * nk + Fy.^2 * nk * nk);

%去除重复的iq，并输出费米面附近iqx，iqy，iqz三个值。保存在vq5当中，可直接用vq5散点图画费米面
vq = [vq2,vq3];
vq5 = [];
vqq = [];
for ip = 1:length(vq2)
    iq = vq2(ip);
    iqx = mod(iq-1,nk) + 1;
    iqy = floor( (mod(iq-1,nk * nk))/nk ) + 1;
    iqz = floor( (iq-1)/(nk * nk) ) + 1;
    iqq = [iqx,iqy,iqz]';
    vq5 = [vq5,iqq];
    iqx = mod((1-iqx+nk),nk) + 1;
    iqy = mod((1-iqy+nk),nk) + 1;
    iqz = mod((1-iqz+nkz),nkz) + 1;
    iq = (iqz - 1) * nk * nk + (iqy-1) * nk + iqx;
    vqq = [vqq,iq];
end
for ip = 1:length(vq3)
    iq = vq3(ip);
    iqx = mod(iq-1,nk) + 1;
    iqy = floor( (mod(iq-1,nk * nk))/nk ) + 1;
    iqz = floor( (iq-1)/(nk * nk) ) + 1;
    iqq = [iqx,iqy,iqz]';
    vq5 = [vq5,iqq];
    iqx = mod((1-iqx+nk),nk) + 1;
    iqy = mod((1-iqy+nk),nk) + 1;
    iqz = mod((1-iqz+nkz),nkz) + 1;
    iq = (iqz - 1) * nk * nk + (iqy-1) * nk + iqx;
    vqq = [vqq,iq];
end
vq5 = vq5';
vq4 = [vq2,vq3];

I = eye(16);
U = zeros(16);
U(1,1) = u;U(6,6) = u;U(11,11) = u;U(16,16) = u;
gams = cell(nq,1);
gamc = cell(nq,1);
for iq = 1:nq               
    chis = inv(I - chi{iq} * U) * chi{iq};
    chic = inv(I + chi{iq} * U) * chi{iq};
    gams{iq} = 0.25 * U * (3 * chis - chic) * U;
    gamc{iq} = 0.25 * U * (3 * chis + chic) * U;
end
chis = zeros(length(vq5));
chic = zeros(length(vq5));
%计算超导能隙矩阵
for ip = 1:length(vq5)
    alpha = 2;
    if ip>length(vq2)
        alpha = 3;
    end
    for iq = ip:length(vq5)
        beta = 2;
        if iq>length(vq2)
            beta = 3;
        end
        iqx = mod((vq5(ip,1) - vq5(iq,1) + nk),nk) + 1;
        iqy = mod((vq5(ip,2) - vq5(iq,2) + nk),nk) + 1;
        iqz = mod((vq5(ip,3) - vq5(iq,3) + nkz),nkz) + 1;
        nq1 = (iqz - 1) * nk * nk + (iqy-1) * nk + iqx;
        
        iqx = mod((vq5(ip,1) + vq5(iq,1) + nk),nk) + 1;
        iqy = mod((vq5(ip,2) + vq5(iq,2) + nk),nk) + 1;
        iqz = mod((vq5(ip,3) + vq5(iq,3) + nkz),nkz) + 1;
        nq2 = (iqz - 1) * nk * nk + (iqy-1) * nk + iqx;  
        
        chi0 = 0;
        chi1 = 0;
        for l1 = 1:4
            for l2 = 1:4
                for l3 = 1:4
                    for l4 = 1:4
                        chi0 = chi0 + (U(l1 * 4 + l2 - 4,l3 * 4+l4 - 4) + gams{nq1}(l1 * 4 + l2 - 4,l3 * 4+l4 - 4) + gams{nq2}(l1 * 4+l4 - 4,l3 * 4 + l2 - 4))...
                             * (conj(mat{l1 + alpha * 4 - 4}(vq4(ip))) * conj(mat{l3 + alpha * 4 - 4}(vqq(ip)))...
                            *mat{l4 + beta * 4 - 4,1}(vqq(iq))*mat{l2 + beta * 4 - 4}(vq4(iq)));
                        
                        chi1 = chi1 + (-gamc{nq1}(l1 * 4 + l2 - 4,l3 * 4+l4 - 4) + gamc{nq2}(l1 * 4+l4 - 4,l3 * 4 + l2 - 4))...
                             * (conj(mat{l1 + alpha * 4 - 4}(vq4(ip))) * conj(mat{l3 + alpha * 4 - 4}(vqq(ip)))...
                            *mat{l4 + beta * 4 - 4}(vqq(iq))*mat{l2 + beta * 4 - 4}(vq4(iq)));
                    end
                end
            end
        end
        if iq > length(vq2)
            chis(ip,iq) = real(chi0)/eng3(vq4(iq));
            chic(ip,iq) = real(chi1)/eng3(vq4(iq));
        else
            chis(ip,iq) = real(chi0)/eng2(vq4(iq));
            chic(ip,iq) = real(chi1)/eng2(vq4(iq));
        end
    end
end
chis = (chis + chis'-diag(diag(chis)));
chic = (chic + chic'-diag(diag(chic)));
chis = chis/(length(vq5));
chic = chic/(length(vq5));
end

