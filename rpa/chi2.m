function AA = chi2(nk,nkz)
    % 计算无相互作用时的极化率
    global eng;
    global mat;
    global mu;
    delta = 1e-5;
    eng1 = cell(4,1);
    mat1 = cell(16,1);
    chi = cell((nk)*(nk)*(nkz),1);
    nq = fix(nk/2 + 0.1);
    nqz = nkz;
    for iq = 1:(nq + 1)*(nq + 1)*(nqz)              
        chiq = zeros(16,16);
        iqx = mod(iq - 1,nq + 1) + 1;
        iqy = floor( (mod(iq - 1,(nq + 1)*(nq + 1)))/(nq + 1) ) + 1;
        iqz = floor( (iq - 1)/((nq + 1)*(nq + 1)) ) + 1;
        for j = 1:4
            eng1{j,1} = circshift(eng{j,1},[ - (iqx - 1), - (iqy - 1), - (iqz - 1)]);   
            enge{j,1} = eng{j,1} < mu;
            enge1{j,1} = eng1{j,1} < mu;
            %enge{j,1} = 1./(exp((eng{j,1} - mu).*B) + 1);      
            %enge1{j,1} = 1./(exp((eng1{j,1} - mu).*B) + 1); 
        end
        for j = 1:16
            mat1{j,1} = circshift(mat{j,1},[ - (iqx - 1), - (iqy - 1), - (iqz - 1)]);
        end    
        for iiii = 1:16
            l1 = floor((iiii - 1)/4 + 1);
            l2 = mod(iiii - 1,4) + 1;
            for jjjj = 1:16
                l4 = floor((jjjj - 1)/4 + 1);
                l3 = mod(jjjj - 1,4) + 1;
                chi0 = 0;
                for alpha = 1:4
                    for beta = 1:4
                        chi0 = chi0 + sum(sum(sum( mat{l4 + alpha * 4 - 4}.*conj(mat{l1 + alpha * 4 - 4})...
                            .*mat1{l2 + beta * 4 - 4}.*conj(mat1{l3 + beta * 4 - 4})...
                            .*( enge1{beta,1} - enge{alpha,1} )./(eng{alpha,1} - eng1{beta,1} + 1i*delta) )));
    %                    chi0 = chi0 + sum(sum(sum( mat{l4 + alpha * 4 - 4}.*conj(mat{l1 + alpha * 4 - 4})...
    %                        .*mat1{l2 + beta * 4 - 4}.*conj(mat1{l3 + beta * 4 - 4})...
    %                        .*( enge1{beta,1} - enge{alpha,1} )./(eng{alpha,1} - eng1{beta,1} + 1i*delta) )));
                    end
                end
                chiq(iiii,jjjj) = chi0;
            end
        end
        chi{(iqz - 1) * nk * nk + (iqy - 1)*nk + iqx} = (chiq)/(nk * nk * nkz);
        if iqx > 1.1
            chi{(iqz - 1) * nk * nk + (iqy - 1)*nk + nk + 2 - iqx} = chi{(iqz - 1) * nk * nk + (iqy - 1)*nk + iqx};
        end
        if iqy > 1.1
            chi{(iqz - 1) * nk * nk + (nk + 1 - iqy)*nk + iqx} = chi{(iqz - 1) * nk * nk + (iqy - 1)*nk + iqx};
        end
        if iqx > 1.1 && iqy>1.1
            chi{(iqz - 1) * nk * nk + (nk + 1 - iqy)*nk + nk + 2 - iqx} = chi{(iqz - 1) * nk * nk + (iqy - 1)*nk + iqx};
        end
    end
    AA = chi;
end


