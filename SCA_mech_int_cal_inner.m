function [mech_int_part] = SCA_mech_int_cal_inner(m, n, p, T, inclusionlist_total, matrixlist, Ni_total, Nm, idx_i_total, idx_m, idx)
% 用于计算multiphase的交互张量
% T为6*6的cell数组
% inclusionlist_total汇总了所有的inclusionlist的数据, Ni_total同样汇总了所有Ni的数据
% 先将cell化的T转化为矩阵
T11 = T{1,1}; T12 = T{1,2}; T13 = T{1,3}; T14 = T{1,4}; T15 = T{1,5}; T16 = T{1,6};
T21 = T{2,1}; T22 = T{2,2}; T23 = T{2,3}; T24 = T{2,4}; T25 = T{2,5}; T26 = T{2,6};
T31 = T{3,1}; T32 = T{3,2}; T33 = T{3,3}; T34 = T{3,4}; T35 = T{3,5}; T36 = T{3,6};
T41 = T{4,1}; T42 = T{4,2}; T43 = T{4,3}; T44 = T{4,4}; T45 = T{4,5}; T46 = T{4,6};
T51 = T{5,1}; T52 = T{5,2}; T53 = T{5,3}; T54 = T{5,4}; T55 = T{5,5}; T56 = T{5,6};
T61 = T{6,1}; T62 = T{6,2}; T63 = T{6,3}; T64 = T{6,4}; T65 = T{6,5}; T66 = T{6,6};
% 还原inclusionlist和Ni的值
Nt = length(Ni_total) + 1; % 计算分簇数目
Nt_total = [Ni_total, Nm]; % 计算所有簇的值

% 定义所有的交互张量，采用cell数组存储
D_total = cell(Nt, Nt);
for i = 1:Nt
    for j = 1:Nt
        D_total{i,j} = zeros(6*Nt_total(i), 6*Nt_total(j));
    end
end

% 依次计算交互张量
for total_Nc = 1:Nt
    if total_Nc ~= Nt
        display(['start calculate interaction tensor in ',num2str(total_Nc),'-th inclusion']);
        Nmain = Ni_total(total_Nc); % 选定当前计算的phase，并得到其分簇数目
        %Nsub = [Ni_total, Nm];
        %Nsub(total_Nc) = []; % 得到其他phase的分簇数目
    else
        display(['start calculate interaction tensor in matrix']);
        Nmain = Nm; % 选定当前计算的phase(即为matrix)，并得到其分簇数目
        %Nsub = [Ni_total, Nm];
        %Nsub = Ni_total; % 得到其他phase的分簇数目(即为所有inclusion)
    end
%     % 构建当前phase与其他phase间的交互张量
%     D_current = zeros(Nt,)
    for i = 1:Nmain
        if total_Nc ~= Nt
            display(['calculating ',num2str(i),'-th clusters in ',num2str(total_Nc),'-th inclusion']);
            %phaselist = inclusionlist_total(:, total_Nc);
            phaselist = inclusionlist_total{total_Nc, 1};
            %idx_phase = idx_i_total{1,total_Nc};
            idx_phase = idx_i_total{total_Nc,1};
        else
            display(['calculating ',num2str(i),'-th clusters in matrix']);
            phaselist = matrixlist;
            idx_phase = idx_m;
        end
        sigmavector = zeros(m*n*p,1);
        sigmavector(phaselist) = (idx_phase==i);
        sigma = reshape(sigmavector,m,n,p);
        fftsigma = fftn(sigma);
        fftsigma = fftshift(fftsigma);
        % FFT transformation
        fftD11 = T11.*fftsigma; fftD12 = T12.*fftsigma; fftD13 = T13.*fftsigma;fftD14 = T14.*fftsigma; fftD15 = T15.*fftsigma; fftD16 = T16.*fftsigma;
        fftD21 = T21.*fftsigma; fftD22 = T22.*fftsigma; fftD23 = T23.*fftsigma;fftD24 = T24.*fftsigma; fftD25 = T25.*fftsigma; fftD26 = T26.*fftsigma;
        fftD31 = T31.*fftsigma; fftD32 = T32.*fftsigma; fftD33 = T33.*fftsigma;fftD34 = T34.*fftsigma; fftD35 = T35.*fftsigma; fftD36 = T36.*fftsigma;
        fftD41 = T41.*fftsigma; fftD42 = T42.*fftsigma; fftD43 = T43.*fftsigma;fftD44 = T44.*fftsigma; fftD45 = T45.*fftsigma; fftD46 = T46.*fftsigma;
        fftD51 = T51.*fftsigma; fftD52 = T52.*fftsigma; fftD53 = T53.*fftsigma;fftD54 = T54.*fftsigma; fftD55 = T55.*fftsigma; fftD56 = T56.*fftsigma;
        fftD61 = T61.*fftsigma; fftD62 = T62.*fftsigma; fftD63 = T63.*fftsigma;fftD64 = T64.*fftsigma; fftD65 = T65.*fftsigma; fftD66 = T66.*fftsigma;
        % Reverse back
        D11 = ifftn(ifftshift(fftD11)); D12 = ifftn(ifftshift(fftD12)); D13 = ifftn(ifftshift(fftD13));D14 = ifftn(ifftshift(fftD14)); D15 = ifftn(ifftshift(fftD15)); D16 = ifftn(ifftshift(fftD16));
        D21 = ifftn(ifftshift(fftD21)); D22 = ifftn(ifftshift(fftD22)); D23 = ifftn(ifftshift(fftD23));D24 = ifftn(ifftshift(fftD24)); D25 = ifftn(ifftshift(fftD25)); D26 = ifftn(ifftshift(fftD26));
        D31 = ifftn(ifftshift(fftD31)); D32 = ifftn(ifftshift(fftD32)); D33 = ifftn(ifftshift(fftD33));D34 = ifftn(ifftshift(fftD34)); D35 = ifftn(ifftshift(fftD35)); D36 = ifftn(ifftshift(fftD36));
        D41 = ifftn(ifftshift(fftD41)); D42 = ifftn(ifftshift(fftD42)); D43 = ifftn(ifftshift(fftD43));D44 = ifftn(ifftshift(fftD44)); D45 = ifftn(ifftshift(fftD45)); D46 = ifftn(ifftshift(fftD46));
        D51 = ifftn(ifftshift(fftD51)); D52 = ifftn(ifftshift(fftD52)); D53 = ifftn(ifftshift(fftD53));D54 = ifftn(ifftshift(fftD54)); D55 = ifftn(ifftshift(fftD55)); D56 = ifftn(ifftshift(fftD56));
        D61 = ifftn(ifftshift(fftD61)); D62 = ifftn(ifftshift(fftD62)); D63 = ifftn(ifftshift(fftD63));D64 = ifftn(ifftshift(fftD64)); D65 = ifftn(ifftshift(fftD65)); D66 = ifftn(ifftshift(fftD66));
        % If m, n, l are even, D_ij (i ~= j) are complex numbers, choose the real part
        D11 = real(reshape(D11,[],1)); D12 = real(reshape(D12,[],1)); D13 = real(reshape(D13,[],1));D14 = real(reshape(D14,[],1)); D15 = real(reshape(D15,[],1)); D16 = real(reshape(D16,[],1));
        D21 = real(reshape(D21,[],1)); D22 = real(reshape(D22,[],1)); D23 = real(reshape(D23,[],1));D24 = real(reshape(D24,[],1)); D25 = real(reshape(D25,[],1)); D26 = real(reshape(D26,[],1));
        D31 = real(reshape(D31,[],1)); D32 = real(reshape(D32,[],1)); D33 = real(reshape(D33,[],1));D34 = real(reshape(D34,[],1)); D35 = real(reshape(D35,[],1)); D36 = real(reshape(D36,[],1));
        D41 = real(reshape(D41,[],1)); D42 = real(reshape(D42,[],1)); D43 = real(reshape(D43,[],1));D44 = real(reshape(D44,[],1)); D45 = real(reshape(D45,[],1)); D46 = real(reshape(D46,[],1));
        D51 = real(reshape(D51,[],1)); D52 = real(reshape(D52,[],1)); D53 = real(reshape(D53,[],1));D54 = real(reshape(D54,[],1)); D55 = real(reshape(D55,[],1)); D56 = real(reshape(D56,[],1));
        D61 = real(reshape(D61,[],1)); D62 = real(reshape(D62,[],1)); D63 = real(reshape(D63,[],1));D64 = real(reshape(D64,[],1)); D65 = real(reshape(D65,[],1)); D66 = real(reshape(D66,[],1));
        
        for Nsub = 1:Nt % 计算交互张量
            if Nsub ~= Nt
                %sub_list = inclusionlist_total(:, Nsub);
                sub_list = inclusionlist_total{Nsub, 1};
            else
                sub_list = matrixlist;
            end
%             temp_D = zeros(6*Nmain, 6*Nt_total(Nsub));
            for j = 1:Nt_total(Nsub)
                Cindex = sub_list&(idx==j);
                Dtemp = [sum(D11(Cindex)),sum(D12(Cindex)),sum(D13(Cindex)),sum(D14(Cindex)),sum(D15(Cindex)),sum(D16(Cindex));...
                         sum(D21(Cindex)),sum(D22(Cindex)),sum(D23(Cindex)),sum(D24(Cindex)),sum(D25(Cindex)),sum(D26(Cindex));...
                         sum(D31(Cindex)),sum(D32(Cindex)),sum(D33(Cindex)),sum(D34(Cindex)),sum(D35(Cindex)),sum(D36(Cindex));...
                         sum(D41(Cindex)),sum(D42(Cindex)),sum(D43(Cindex)),sum(D44(Cindex)),sum(D45(Cindex)),sum(D46(Cindex));...
                         sum(D51(Cindex)),sum(D52(Cindex)),sum(D53(Cindex)),sum(D54(Cindex)),sum(D55(Cindex)),sum(D56(Cindex));...
                         sum(D61(Cindex)),sum(D62(Cindex)),sum(D63(Cindex)),sum(D64(Cindex)),sum(D65(Cindex)),sum(D66(Cindex))];
                %temp_D(6*j-5:6*j,6*i-5:6*i) = Dtemp/sum(Cindex);
                D_total{Nsub, total_Nc}(6*j-5:6*j,6*i-5:6*i) = Dtemp/sum(Cindex);
            end
        end     
    end
end
mech_int_part = cell2mat(D_total);
end
        %display([num2str(i),'-th clusters in inclusion']);
    