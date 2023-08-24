function [mech_int_part1, mech_int_part2] = SCA_mech_int_cal(m, n, p, em, en, ep, inclusionlist_total, matrixlist, Ni_total, Nm, idx_i_total, idx_m, idx)
% 计算交互张量的两个部分part1和part2
%% 计算第一部分
% First part of the interaction tensor ---------------------------------------
display('Start computing the first part of the mech interaction matrice');
% 更新Green函数的第一部分
T11 = zeros(m,n,p);T12 = zeros(m,n,p);T13 = zeros(m,n,p);T14 = zeros(m,n,p);T15 = zeros(m,n,p);T16 = zeros(m,n,p);
                   T22 = zeros(m,n,p);T23 = zeros(m,n,p);T24 = zeros(m,n,p);T25 = zeros(m,n,p);T26 = zeros(m,n,p);
                                      T33 = zeros(m,n,p);T34 = zeros(m,n,p);T35 = zeros(m,n,p);T36 = zeros(m,n,p);
                                                         T44 = zeros(m,n,p);T45 = zeros(m,n,p);T46 = zeros(m,n,p);
                                                                            T55 = zeros(m,n,p);T56 = zeros(m,n,p);
                                                                                               T66 = zeros(m,n,p);
for i = 1:m
    for j = 1:n
        for k = 1:p
            e1 = em(i); e2 = en(j); e3 = ep(k); E = norm([e1, e2, e3]);
            T11(i,j,k) = 1/E^2*(4*e1^2);
            T12(i,j,k) = 0;
            T13(i,j,k) = 0;
            T14(i,j,k) = 1/E^2*(2*e1*e2);
            T15(i,j,k) = 0;
            T16(i,j,k) = 1/E^2*(2*e1*e3);
            T22(i,j,k) = 1/E^2*(4*e2^2);
            T23(i,j,k) = 0;
            T24(i,j,k) = 1/E^2*(2*e1*e2);
            T25(i,j,k) = 1/E^2*(2*e2*e3);
            T26(i,j,k) = 0;
            T33(i,j,k) = 1/E^2*(4*e3^2);
            T34(i,j,k) = 0;
            T35(i,j,k) = 1/E^2*(2*e2*e3);
            T36(i,j,k) = 1/E^2*(2*e1*e3);
            T44(i,j,k) = 1/E^2*(e1^2+e2^2);
            T45(i,j,k) = 1/E^2*(e1*e3);
            T46(i,j,k) = 1/E^2*(e2*e3);
            T55(i,j,k) = 1/E^2*(e2^2+e3^2);
            T56(i,j,k) = 1/E^2*(e1*e2);
            T66(i,j,k) = 1/E^2*(e1^2+e3^2);
        end
    end
end
% 修改Green函数
T14 = T14*2; T15 = T15*2; T16 = T16*2; T24 = T24*2; T25 = T25*2; T26 = T26*2;
T34 = T34*2; T35 = T35*2; T36 = T36*2; 
T44 = T44*4; T45 = T45*4; T46 = T46*4; T55 = T55*4; T56 = T56*4; T66 = T66*4;
% 对称
T21 = T12; T31 = T13; T41 = T14; T51 = T15; T61 = T16;
T32 = T23; T42 = T24; T52 = T25; T62 = T26;
T43 = T34; T53 = T35; T63 = T36;
T54 = T45; T64 = T46;
T65 = T56;
% 修改Green函数在0频率处的值
cpx = ceil((m+1)/2); cpy = ceil((n+1)/2); cpz = ceil((p+1)/2);
T11(cpx,cpy,cpz) = 0; T12(cpx,cpy,cpz) = 0; T13(cpx,cpy,cpz) = 0; T14(cpx,cpy,cpz) = 0; T15(cpx,cpy,cpz) = 0; T16(cpx,cpy,cpz) = 0;
T21(cpx,cpy,cpz) = 0; T22(cpx,cpy,cpz) = 0; T23(cpx,cpy,cpz) = 0; T24(cpx,cpy,cpz) = 0; T25(cpx,cpy,cpz) = 0; T26(cpx,cpy,cpz) = 0;
T31(cpx,cpy,cpz) = 0; T32(cpx,cpy,cpz) = 0; T33(cpx,cpy,cpz) = 0; T34(cpx,cpy,cpz) = 0; T35(cpx,cpy,cpz) = 0; T36(cpx,cpy,cpz) = 0;
T41(cpx,cpy,cpz) = 0; T42(cpx,cpy,cpz) = 0; T43(cpx,cpy,cpz) = 0; T44(cpx,cpy,cpz) = 0; T45(cpx,cpy,cpz) = 0; T46(cpx,cpy,cpz) = 0;
T51(cpx,cpy,cpz) = 0; T52(cpx,cpy,cpz) = 0; T53(cpx,cpy,cpz) = 0; T54(cpx,cpy,cpz) = 0; T55(cpx,cpy,cpz) = 0; T56(cpx,cpy,cpz) = 0;
T61(cpx,cpy,cpz) = 0; T62(cpx,cpy,cpz) = 0; T63(cpx,cpy,cpz) = 0; T64(cpx,cpy,cpz) = 0; T65(cpx,cpy,cpz) = 0; T66(cpx,cpy,cpz) = 0;

% Define all the interaction tensors
% 定义所有的交互张量？
% 将交互张量的计算写成函数计算，5相材料
% 需要将Green函数进行cell化
T_total = cell(6,6);
T_total{1,1} = T11; T_total{1,2} = T12; T_total{1,3} = T13; T_total{1,4} = T14; T_total{1,5} = T15; T_total{1,6} = T16; 
T_total{2,1} = T21; T_total{2,2} = T22; T_total{2,3} = T23; T_total{2,4} = T24; T_total{2,5} = T25; T_total{2,6} = T26; 
T_total{3,1} = T31; T_total{3,2} = T32; T_total{3,3} = T33; T_total{3,4} = T34; T_total{3,5} = T35; T_total{3,6} = T36; 
T_total{4,1} = T41; T_total{4,2} = T42; T_total{4,3} = T43; T_total{4,4} = T44; T_total{4,5} = T45; T_total{4,6} = T46; 
T_total{5,1} = T51; T_total{5,2} = T52; T_total{5,3} = T53; T_total{5,4} = T54; T_total{5,5} = T55; T_total{5,6} = T56; 
T_total{6,1} = T61; T_total{6,2} = T62; T_total{6,3} = T63; T_total{6,4} = T64; T_total{6,5} = T65; T_total{6,6} = T66; 
% % 之后将不同inclusionlist写在一起
% inclusionlist_total = [inclusionlist1, inclusionlist2, inclusionlist3, inclusionlist4];
% Ni_total = [Ni1,Ni2,Ni3,Ni4];
[mech_int_part1] = SCA_mech_int_cal_inner(m, n, p, T_total, inclusionlist_total, matrixlist, Ni_total, Nm, idx_i_total, idx_m, idx);

%% 计算第二部分
% Second part of the interaction tensor ---------------------------------------
display('Start computing the second part of the mech interaction matrice');
T11 = zeros(m,n,p);T12 = zeros(m,n,p);T13 = zeros(m,n,p);T14 = zeros(m,n,p);T15 = zeros(m,n,p);T16 = zeros(m,n,p);
                   T22 = zeros(m,n,p);T23 = zeros(m,n,p);T24 = zeros(m,n,p);T25 = zeros(m,n,p);T26 = zeros(m,n,p);
                                      T33 = zeros(m,n,p);T34 = zeros(m,n,p);T35 = zeros(m,n,p);T36 = zeros(m,n,p);
                                                         T44 = zeros(m,n,p);T45 = zeros(m,n,p);T46 = zeros(m,n,p);
                                                                            T55 = zeros(m,n,p);T56 = zeros(m,n,p);
                                                                                               T66 = zeros(m,n,p);
% Green函数右边部分需要修改，7.31
for i = 1:m
    for j = 1:n
        for k = 1:p
            e1 = em(i); e2 = en(j); e3 = ep(k); E = norm([e1, e2, e3]);
            T11(i,j,k) = -e1^4/E^4;
            T12(i,j,k) = -e1^2*e2^2/E^4;
            T13(i,j,k) = -e1^2*e3^2/E^4;
            T14(i,j,k) = -e1^3*e2/E^4;
            T15(i,j,k) = -e1^2*e2*e3/E^4;
            T16(i,j,k) = -e1^3*e3/E^4;
            T22(i,j,k) = -e2^4/E^4;
            T23(i,j,k) = -e2^2*e3^2/E^4;
            T24(i,j,k) = -e2^3*e1/E^4;
            T25(i,j,k) = -e2^3*e3/E^4;
            T26(i,j,k) = -e2^2*e1*e3/E^4;
            T33(i,j,k) = -e3^4/E^4;
            T34(i,j,k) = -e3^2*e1*e2/E^4;
            T35(i,j,k) = -e3^3*e2/E^4;
            T36(i,j,k) = -e3^3*e1/E^4;
            T44(i,j,k) = -e1^2*e2^2/E^4;
            T45(i,j,k) = -e1*e2^2*e3/E^4;
            T46(i,j,k) = -e1^2*e2*e3/E^4;
            T55(i,j,k) = -e2^2*e3^2/E^4;
            T56(i,j,k) = -e1*e2*e3^2/E^4;
            T66(i,j,k) = -e1^2*e3^2/E^4;
        end
    end
end
% 修改Green函数
T14 = T14*2; T15 = T15*2; T16 = T16*2; T24 = T24*2; T25 = T25*2; T26 = T26*2;
T34 = T34*2; T35 = T35*2; T36 = T36*2; 
T44 = T44*4; T45 = T45*4; T46 = T46*4; T55 = T55*4; T56 = T56*4; T66 = T66*4;
% 对称
T21 = T12; T31 = T13; T41 = T14; T51 = T15; T61 = T16;
T32 = T23; T42 = T24; T52 = T25; T62 = T26;
T43 = T34; T53 = T35; T63 = T36;
T54 = T45; T64 = T46;
T65 = T56;
% 修改Green函数在0频率处的值
cpx = ceil((m+1)/2); cpy = ceil((n+1)/2); cpz = ceil((p+1)/2);
T11(cpx,cpy,cpz) = 0; T12(cpx,cpy,cpz) = 0; T13(cpx,cpy,cpz) = 0; T14(cpx,cpy,cpz) = 0; T15(cpx,cpy,cpz) = 0; T16(cpx,cpy,cpz) = 0;
T21(cpx,cpy,cpz) = 0; T22(cpx,cpy,cpz) = 0; T23(cpx,cpy,cpz) = 0; T24(cpx,cpy,cpz) = 0; T25(cpx,cpy,cpz) = 0; T26(cpx,cpy,cpz) = 0;
T31(cpx,cpy,cpz) = 0; T32(cpx,cpy,cpz) = 0; T33(cpx,cpy,cpz) = 0; T34(cpx,cpy,cpz) = 0; T35(cpx,cpy,cpz) = 0; T36(cpx,cpy,cpz) = 0;
T41(cpx,cpy,cpz) = 0; T42(cpx,cpy,cpz) = 0; T43(cpx,cpy,cpz) = 0; T44(cpx,cpy,cpz) = 0; T45(cpx,cpy,cpz) = 0; T46(cpx,cpy,cpz) = 0;
T51(cpx,cpy,cpz) = 0; T52(cpx,cpy,cpz) = 0; T53(cpx,cpy,cpz) = 0; T54(cpx,cpy,cpz) = 0; T55(cpx,cpy,cpz) = 0; T56(cpx,cpy,cpz) = 0;
T61(cpx,cpy,cpz) = 0; T62(cpx,cpy,cpz) = 0; T63(cpx,cpy,cpz) = 0; T64(cpx,cpy,cpz) = 0; T65(cpx,cpy,cpz) = 0; T66(cpx,cpy,cpz) = 0;

% 重定义T_total
T_total = cell(6,6);
T_total{1,1} = T11; T_total{1,2} = T12; T_total{1,3} = T13; T_total{1,4} = T14; T_total{1,5} = T15; T_total{1,6} = T16; 
T_total{2,1} = T21; T_total{2,2} = T22; T_total{2,3} = T23; T_total{2,4} = T24; T_total{2,5} = T25; T_total{2,6} = T26; 
T_total{3,1} = T31; T_total{3,2} = T32; T_total{3,3} = T33; T_total{3,4} = T34; T_total{3,5} = T35; T_total{3,6} = T36; 
T_total{4,1} = T41; T_total{4,2} = T42; T_total{4,3} = T43; T_total{4,4} = T44; T_total{4,5} = T45; T_total{4,6} = T46; 
T_total{5,1} = T51; T_total{5,2} = T52; T_total{5,3} = T53; T_total{5,4} = T54; T_total{5,5} = T55; T_total{5,6} = T56; 
T_total{6,1} = T61; T_total{6,2} = T62; T_total{6,3} = T63; T_total{6,4} = T64; T_total{6,5} = T65; T_total{6,6} = T66; 
[mech_int_part2] = SCA_mech_int_cal_inner(m, n, p, T_total, inclusionlist_total, matrixlist, Ni_total, Nm, idx_i_total, idx_m, idx);

end
