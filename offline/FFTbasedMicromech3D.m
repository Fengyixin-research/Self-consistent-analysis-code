function [etotal_11,etotal_22,etotal_33,etotal_12,etotal_23,etotal_13] = FFTbasedMicromech3D(phase, material, C0, dE_initial)
% Input
% phase: material phase, 0/1/2/3/4 for five material phases in composite
% material:lamda,mu
% material: [lamda, mu, material_property], first inclusion, then matrix
% C0: prescribed homogenized stiffness tensor
% output
% the strian at different point, direction is 11,22,33,12,23,13 in 3D

% number of material phases
N_phase = size(material(:,end),1);
C_total = cell(N_phase, 1);
for i = 1:N_phase
    C_total{i,1} = C_cal_3D_orthotropic_comsol(material(i,1), material(i,2), material(i,3), material(i,4), material(i,5), material(i,6), material(i,7), material(i,8), material(i,9));
end

% reference material
lamda_0 = C0(1,2); mu_0 = C0(4,4); 
coef_1 = 1/(4*mu_0); coef_2 = (lamda_0+mu_0)/(mu_0*(lamda_0+2*mu_0));
% 3D green function,6*6*m*n*p
[m,n,p] = size(phase);
% frequency
em = (ceil(-m/2):floor((m-1)/2));
en = (ceil(-n/2):floor((n-1)/2));
ep = (ceil(-p/2):floor((p-1)/2));
G = Green_cal(coef_1, coef_2, em, en, ep, m, n, p);
G_11 = G{1,1}; G_12 = G{1,2}; G_13 = G{1,3}; G_14 = G{1,4}; G_15 = G{1,5}; G_16 = G{1,6};
G_21 = G{2,1}; G_22 = G{2,2}; G_23 = G{2,3}; G_24 = G{2,4}; G_25 = G{2,5}; G_26 = G{2,6};
G_31 = G{3,1}; G_32 = G{3,2}; G_33 = G{3,3}; G_34 = G{3,4}; G_35 = G{3,5}; G_36 = G{3,6};
G_41 = G{4,1}; G_42 = G{4,2}; G_43 = G{4,3}; G_44 = G{4,4}; G_45 = G{4,5}; G_46 = G{4,6};
G_51 = G{5,1}; G_52 = G{5,2}; G_53 = G{5,3}; G_54 = G{5,4}; G_55 = G{5,5}; G_56 = G{5,6};
G_61 = G{6,1}; G_62 = G{6,2}; G_63 = G{6,3}; G_64 = G{6,4}; G_65 = G{6,5}; G_66 = G{6,6};

%% initialization with the zero values
etotal_11=zeros(m,n,p);etotal_22=zeros(m,n,p);etotal_33=zeros(m,n,p);etotal_12=zeros(m,n,p);etotal_23=zeros(m,n,p);etotal_13=zeros(m,n,p);
sigma_11=zeros(m,n,p);sigma_22=zeros(m,n,p);sigma_33=zeros(m,n,p);sigma_12=zeros(m,n,p);sigma_23=zeros(m,n,p);sigma_13=zeros(m,n,p);
% polar stresses
sigma_polar_11=zeros(m,n,p);sigma_polar_22=zeros(m,n,p);sigma_polar_33=zeros(m,n,p);sigma_polar_12=zeros(m,n,p);sigma_polar_23=zeros(m,n,p);sigma_polar_13=zeros(m,n,p);

disp(['Start calculate FFT_mech preprocess.']);
% incremental of total strain
dE = zeros(6,1);
dE_total = dE + dE_initial;
etotal_11 = etotal_11 + dE_total(1); etotal_22 = etotal_22 + dE_total(2); etotal_33 = etotal_33 + dE_total(3);
etotal_12 = etotal_12 + dE_total(4); etotal_23 = etotal_23 + dE_total(5); etotal_13 = etotal_13 + dE_total(6);
iteration = 0;
while  1
    etotal_11_initial = etotal_11; etotal_22_initial = etotal_22; etotal_33_initial = etotal_33;
    etotal_12_initial = etotal_12; etotal_23_initial = etotal_23; etotal_13_initial = etotal_13; 
    for i = 1:m
        for j = 1:n
            for k = 1:p
                temp_C = C_total{phase(i,j,k)+1,1};
                temp_sigma = (temp_C - C0)*[etotal_11(i,j,k);etotal_22(i,j,k);etotal_33(i,j,k);etotal_12(i,j,k);etotal_23(i,j,k);etotal_13(i,j,k)];
                sigma_polar_11(i,j,k) = temp_sigma(1);
                sigma_polar_22(i,j,k) = temp_sigma(2);
                sigma_polar_33(i,j,k) = temp_sigma(3);
                sigma_polar_12(i,j,k) = temp_sigma(4);
                sigma_polar_23(i,j,k) = temp_sigma(5);
                sigma_polar_13(i,j,k) = temp_sigma(6);
            end
        end
    end
     % FFT transform
    sigma_polarfft_11 = fftshift(fftn(sigma_polar_11));sigma_polarfft_22 = fftshift(fftn(sigma_polar_22));sigma_polarfft_33 = fftshift(fftn(sigma_polar_33));
    sigma_polarfft_12 = fftshift(fftn(sigma_polar_12));sigma_polarfft_23 = fftshift(fftn(sigma_polar_23));sigma_polarfft_13 = fftshift(fftn(sigma_polar_13));
    % update strain
    etotalfft_11 = -(G_11.*sigma_polarfft_11 + G_12.*sigma_polarfft_22 + G_13.*sigma_polarfft_33 + G_14.*sigma_polarfft_12 + G_15.*sigma_polarfft_23 + G_16.*sigma_polarfft_13);
    etotalfft_22 = -(G_21.*sigma_polarfft_11 + G_22.*sigma_polarfft_22 + G_23.*sigma_polarfft_33 + G_24.*sigma_polarfft_12 + G_25.*sigma_polarfft_23 + G_26.*sigma_polarfft_13);
    etotalfft_33 = -(G_31.*sigma_polarfft_11 + G_32.*sigma_polarfft_22 + G_33.*sigma_polarfft_33 + G_34.*sigma_polarfft_12 + G_35.*sigma_polarfft_23 + G_36.*sigma_polarfft_13);
    etotalfft_12 = -(G_41.*sigma_polarfft_11 + G_42.*sigma_polarfft_22 + G_43.*sigma_polarfft_33 + G_44.*sigma_polarfft_12 + G_45.*sigma_polarfft_23 + G_46.*sigma_polarfft_13);
    etotalfft_23 = -(G_51.*sigma_polarfft_11 + G_52.*sigma_polarfft_22 + G_53.*sigma_polarfft_33 + G_54.*sigma_polarfft_12 + G_55.*sigma_polarfft_23 + G_56.*sigma_polarfft_13);
    etotalfft_13 = -(G_61.*sigma_polarfft_11 + G_62.*sigma_polarfft_22 + G_63.*sigma_polarfft_33 + G_64.*sigma_polarfft_12 + G_65.*sigma_polarfft_23 + G_66.*sigma_polarfft_13); 
    % inverse FFT
    etotal_11 = real(ifftn(ifftshift(etotalfft_11))); etotal_22 = real(ifftn(ifftshift(etotalfft_22))); etotal_33 = real(ifftn(ifftshift(etotalfft_33)));
    etotal_12 = real(ifftn(ifftshift(etotalfft_12))); etotal_23 = real(ifftn(ifftshift(etotalfft_23))); etotal_13 = real(ifftn(ifftshift(etotalfft_13)));
    etotal_11 = dE_total(1) + etotal_11; etotal_22 = dE_total(2) + etotal_22; etotal_33 = dE_total(3) + etotal_33;
    etotal_12 = dE_total(4) + etotal_12; etotal_23 = dE_total(5) + etotal_23; etotal_13 = dE_total(6) + etotal_13;         
    % update stress
    for i = 1:m
        for j = 1:n
            for k = 1:p
                temp_C = C_total{phase(i,j,k)+1,1};
                temp_sigma = temp_C*[etotal_11(i,j,k);etotal_22(i,j,k);etotal_33(i,j,k);...
                    etotal_12(i,j,k); etotal_23(i,j,k); etotal_13(i,j,k)];
                sigma_11(i,j,k) = temp_sigma(1);
                sigma_22(i,j,k) = temp_sigma(2);
                sigma_33(i,j,k) = temp_sigma(3);
                sigma_12(i,j,k) = temp_sigma(4);
                sigma_23(i,j,k) = temp_sigma(5);
                sigma_13(i,j,k) = temp_sigma(6);
            end
        end
    end              
    er = (sum(sqrt((etotal_11-etotal_11_initial).^2 + (etotal_22-etotal_22_initial).^2 + (etotal_33-etotal_33_initial).^2 + (etotal_12-etotal_12_initial).^2 + ...
        (etotal_23-etotal_23_initial).^2 + (etotal_13-etotal_13_initial).^2), 'all')/m/n/p) / sum(sqrt(dE_total.^2),'all');
    if er < 1e-10
        break;
    end
    iteration = iteration + 1;    
end
end

function [G] = Green_cal(c1, c2, em, en, ep, m, n, p)
% construct the Green function
% first part
Gl_11 = zeros(m,n,p);Gl_12 = zeros(m,n,p);Gl_13 = zeros(m,n,p);Gl_14 = zeros(m,n,p);Gl_15 = zeros(m,n,p);Gl_16 = zeros(m,n,p);
                     Gl_22 = zeros(m,n,p);Gl_23 = zeros(m,n,p);Gl_24 = zeros(m,n,p);Gl_25 = zeros(m,n,p);Gl_26 = zeros(m,n,p);
                                          Gl_33 = zeros(m,n,p);Gl_34 = zeros(m,n,p);Gl_35 = zeros(m,n,p);Gl_36 = zeros(m,n,p);
                                                               Gl_44 = zeros(m,n,p);Gl_45 = zeros(m,n,p);Gl_46 = zeros(m,n,p);
                                                                                    Gl_55 = zeros(m,n,p);Gl_56 = zeros(m,n,p);
                                                                                                         Gl_66 = zeros(m,n,p);
for i = 1:m
    for j = 1:n
        for k = 1:p
            e1 = em(i); e2 = en(j); e3 = ep(k); E = norm([e1, e2, e3]);
            Gl_11(i,j,k) = 1/E^2*(4*e1^2);
            Gl_12(i,j,k) = 0;
            Gl_13(i,j,k) = 0;
            Gl_14(i,j,k) = 1/E^2*(2*e1*e2);
            Gl_15(i,j,k) = 0;
            Gl_16(i,j,k) = 1/E^2*(2*e1*e3);
            Gl_22(i,j,k) = 1/E^2*(4*e2^2);
            Gl_23(i,j,k) = 0;
            Gl_24(i,j,k) = 1/E^2*(2*e1*e2);
            Gl_25(i,j,k) = 1/E^2*(2*e2*e3);
            Gl_26(i,j,k) = 0;
            Gl_33(i,j,k) = 1/E^2*(4*e3^2);
            Gl_34(i,j,k) = 0;
            Gl_35(i,j,k) = 1/E^2*(2*e2*e3);
            Gl_36(i,j,k) = 1/E^2*(2*e1*e3);
            Gl_44(i,j,k) = 1/E^2*(e1^2+e2^2);
            Gl_45(i,j,k) = 1/E^2*(e1*e3);
            Gl_46(i,j,k) = 1/E^2*(e2*e3);
            Gl_55(i,j,k) = 1/E^2*(e2^2+e3^2);
            Gl_56(i,j,k) = 1/E^2*(e1*e2);
            Gl_66(i,j,k) = 1/E^2*(e1^2+e3^2);
        end
    end
end
% Voigt
Gl_14 = Gl_14*2; Gl_15 = Gl_15*2; Gl_16 = Gl_16*2; Gl_24 = Gl_24*2; Gl_25 = Gl_25*2; Gl_26 = Gl_26*2;
Gl_34 = Gl_34*2; Gl_35 = Gl_35*2; Gl_36 = Gl_36*2;
Gl_44 = Gl_44*4; Gl_45 = Gl_45*4; Gl_46 = Gl_46*4; Gl_55 = Gl_55*4; Gl_56 = Gl_56*4; Gl_66 = Gl_66*4;
% symmetry
Gl_21 = Gl_12; 
Gl_31 = Gl_13; Gl_32 = Gl_23;
Gl_41 = Gl_14; Gl_42 = Gl_24; Gl_43 = Gl_34;
Gl_51 = Gl_15; Gl_52 = Gl_25; Gl_53 = Gl_35; Gl_54 = Gl_45;
Gl_61 = Gl_16; Gl_62 = Gl_26; Gl_63 = Gl_36; Gl_64 = Gl_46; Gl_65 = Gl_56;
% zero frequency
cpx = ceil((m+1)/2); cpy = ceil((n+1)/2); cpz = ceil((p+1)/2);
Gl_11(cpx,cpy,cpz) = 0; Gl_12(cpx,cpy,cpz) = 0; Gl_13(cpx,cpy,cpz) = 0; Gl_14(cpx,cpy,cpz) = 0; Gl_15(cpx,cpy,cpz) = 0; Gl_16(cpx,cpy,cpz) = 0;
Gl_21(cpx,cpy,cpz) = 0; Gl_22(cpx,cpy,cpz) = 0; Gl_23(cpx,cpy,cpz) = 0; Gl_24(cpx,cpy,cpz) = 0; Gl_25(cpx,cpy,cpz) = 0; Gl_26(cpx,cpy,cpz) = 0;
Gl_31(cpx,cpy,cpz) = 0; Gl_32(cpx,cpy,cpz) = 0; Gl_33(cpx,cpy,cpz) = 0; Gl_34(cpx,cpy,cpz) = 0; Gl_35(cpx,cpy,cpz) = 0; Gl_36(cpx,cpy,cpz) = 0;
Gl_41(cpx,cpy,cpz) = 0; Gl_42(cpx,cpy,cpz) = 0; Gl_43(cpx,cpy,cpz) = 0; Gl_44(cpx,cpy,cpz) = 0; Gl_45(cpx,cpy,cpz) = 0; Gl_46(cpx,cpy,cpz) = 0; 
Gl_51(cpx,cpy,cpz) = 0; Gl_52(cpx,cpy,cpz) = 0; Gl_53(cpx,cpy,cpz) = 0; Gl_54(cpx,cpy,cpz) = 0; Gl_55(cpx,cpy,cpz) = 0; Gl_56(cpx,cpy,cpz) = 0;
Gl_61(cpx,cpy,cpz) = 0; Gl_62(cpx,cpy,cpz) = 0; Gl_63(cpx,cpy,cpz) = 0; Gl_64(cpx,cpy,cpz) = 0; Gl_65(cpx,cpy,cpz) = 0; Gl_66(cpx,cpy,cpz) = 0;

% second part
Gr_11 = zeros(m,n,p);Gr_12 = zeros(m,n,p);Gr_13 = zeros(m,n,p);Gr_14 = zeros(m,n,p);Gr_15 = zeros(m,n,p);Gr_16 = zeros(m,n,p);
                     Gr_22 = zeros(m,n,p);Gr_23 = zeros(m,n,p);Gr_24 = zeros(m,n,p);Gr_25 = zeros(m,n,p);Gr_26 = zeros(m,n,p);
                                          Gr_33 = zeros(m,n,p);Gr_34 = zeros(m,n,p);Gr_35 = zeros(m,n,p);Gr_36 = zeros(m,n,p);
                                                               Gr_44 = zeros(m,n,p);Gr_45 = zeros(m,n,p);Gr_46 = zeros(m,n,p);
                                                                                    Gr_55 = zeros(m,n,p);Gr_56 = zeros(m,n,p);
                                                                                                         Gr_66 = zeros(m,n,p);
for i = 1:m
    for j = 1:n
        for k = 1:p
            e1 = em(i); e2 = en(j); e3 = ep(k); E = norm([e1, e2, e3]);
            Gr_11(i,j,k) = -e1^4/E^4;
            Gr_12(i,j,k) = -e1^2*e2^2/E^4;
            Gr_13(i,j,k) = -e1^2*e3^2/E^4;
            Gr_14(i,j,k) = -e1^3*e2/E^4;
            Gr_15(i,j,k) = -e1^2*e2*e3/E^4;
            Gr_16(i,j,k) = -e1^3*e3/E^4;
            Gr_22(i,j,k) = -e2^4/E^4;
            Gr_23(i,j,k) = -e2^2*e3^2/E^4;
            Gr_24(i,j,k) = -e2^3*e1/E^4;
            Gr_25(i,j,k) = -e2^3*e3/E^4;
            Gr_26(i,j,k) = -e2^2*e1*e3/E^4;
            Gr_33(i,j,k) = -e3^4/E^4;
            Gr_34(i,j,k) = -e3^2*e1*e2/E^4;
            Gr_35(i,j,k) = -e3^3*e2/E^4;
            Gr_36(i,j,k) = -e3^3*e1/E^4;
            Gr_44(i,j,k) = -e1^2*e2^2/E^4;
            Gr_45(i,j,k) = -e1*e2^2*e3/E^4;
            Gr_46(i,j,k) = -e1^2*e2*e3/E^4;
            Gr_55(i,j,k) = -e2^2*e3^2/E^4;
            Gr_56(i,j,k) = -e1*e2*e3^2/E^4;
            Gr_66(i,j,k) = -e1^2*e3^2/E^4;
        end
    end
end
Gr_14 = Gr_14*2; Gr_15 = Gr_15*2; Gr_16 = Gr_16*2; Gr_24 = Gr_24*2; Gr_25 = Gr_25*2; Gr_26 = Gr_26*2;
Gr_34 = Gr_34*2; Gr_35 = Gr_35*2; Gr_36 = Gr_36*2;
Gr_44 = Gr_44*4; Gr_45 = Gr_45*4; Gr_46 = Gr_46*4; Gr_55 = Gr_55*4; Gr_56 = Gr_56*4; Gr_66 = Gr_66*4;
Gr_21 = Gr_12; 
Gr_31 = Gr_13; Gr_32 = Gr_23;
Gr_41 = Gr_14; Gr_42 = Gr_24; Gr_43 = Gr_34;
Gr_51 = Gr_15; Gr_52 = Gr_25; Gr_53 = Gr_35; Gr_54 = Gr_45;
Gr_61 = Gr_16; Gr_62 = Gr_26; Gr_63 = Gr_36; Gr_64 = Gr_46; Gr_65 = Gr_56;                                                                                                         
Gr_11(cpx,cpy,cpz) = 0; Gr_12(cpx,cpy,cpz) = 0; Gr_13(cpx,cpy,cpz) = 0; Gr_14(cpx,cpy,cpz) = 0; Gr_15(cpx,cpy,cpz) = 0; Gr_16(cpx,cpy,cpz) = 0;
Gr_21(cpx,cpy,cpz) = 0; Gr_22(cpx,cpy,cpz) = 0; Gr_23(cpx,cpy,cpz) = 0; Gr_24(cpx,cpy,cpz) = 0; Gr_25(cpx,cpy,cpz) = 0; Gr_26(cpx,cpy,cpz) = 0;
Gr_31(cpx,cpy,cpz) = 0; Gr_32(cpx,cpy,cpz) = 0; Gr_33(cpx,cpy,cpz) = 0; Gr_34(cpx,cpy,cpz) = 0; Gr_35(cpx,cpy,cpz) = 0; Gr_36(cpx,cpy,cpz) = 0;
Gr_41(cpx,cpy,cpz) = 0; Gr_42(cpx,cpy,cpz) = 0; Gr_43(cpx,cpy,cpz) = 0; Gr_44(cpx,cpy,cpz) = 0; Gr_45(cpx,cpy,cpz) = 0; Gr_46(cpx,cpy,cpz) = 0; 
Gr_51(cpx,cpy,cpz) = 0; Gr_52(cpx,cpy,cpz) = 0; Gr_53(cpx,cpy,cpz) = 0; Gr_54(cpx,cpy,cpz) = 0; Gr_55(cpx,cpy,cpz) = 0; Gr_56(cpx,cpy,cpz) = 0;
Gr_61(cpx,cpy,cpz) = 0; Gr_62(cpx,cpy,cpz) = 0; Gr_63(cpx,cpy,cpz) = 0; Gr_64(cpx,cpy,cpz) = 0; Gr_65(cpx,cpy,cpz) = 0; Gr_66(cpx,cpy,cpz) = 0;


G_11=c1*Gl_11+c2*Gr_11; G_12=c1*Gl_12+c2*Gr_12; G_13=c1*Gl_13+c2*Gr_13; G_14=c1*Gl_14+c2*Gr_14; G_15=c1*Gl_15+c2*Gr_15; G_16=c1*Gl_16+c2*Gr_16;
G_21=c1*Gl_21+c2*Gr_21; G_22=c1*Gl_22+c2*Gr_22; G_23=c1*Gl_23+c2*Gr_23; G_24=c1*Gl_24+c2*Gr_24; G_25=c1*Gl_25+c2*Gr_25; G_26=c1*Gl_26+c2*Gr_26;
G_31=c1*Gl_31+c2*Gr_31; G_32=c1*Gl_32+c2*Gr_32; G_33=c1*Gl_33+c2*Gr_33; G_34=c1*Gl_34+c2*Gr_34; G_35=c1*Gl_35+c2*Gr_35; G_36=c1*Gl_36+c2*Gr_36;
G_41=c1*Gl_41+c2*Gr_41; G_42=c1*Gl_42+c2*Gr_42; G_43=c1*Gl_43+c2*Gr_43; G_44=c1*Gl_44+c2*Gr_44; G_45=c1*Gl_45+c2*Gr_45; G_46=c1*Gl_46+c2*Gr_46;
G_51=c1*Gl_51+c2*Gr_51; G_52=c1*Gl_52+c2*Gr_52; G_53=c1*Gl_53+c2*Gr_53; G_54=c1*Gl_54+c2*Gr_54; G_55=c1*Gl_55+c2*Gr_55; G_56=c1*Gl_56+c2*Gr_56;
G_61=c1*Gl_61+c2*Gr_61; G_62=c1*Gl_62+c2*Gr_62; G_63=c1*Gl_63+c2*Gr_63; G_64=c1*Gl_64+c2*Gr_64; G_65=c1*Gl_65+c2*Gr_65; G_66=c1*Gl_66+c2*Gr_66;

G = cell(6,6);
G{1,1} = G_11; G{1,2} = G_12; G{1,3} = G_13; G{1,4} = G_14; G{1,5} = G_15; G{1,6} = G_16;
G{2,1} = G_21; G{2,2} = G_22; G{2,3} = G_23; G{2,4} = G_24; G{2,5} = G_25; G{2,6} = G_26;
G{3,1} = G_31; G{3,2} = G_32; G{3,3} = G_33; G{3,4} = G_34; G{3,5} = G_35; G{3,6} = G_36;
G{4,1} = G_41; G{4,2} = G_42; G{4,3} = G_43; G{4,4} = G_44; G{4,5} = G_45; G{4,6} = G_46;
G{5,1} = G_51; G{5,2} = G_52; G{5,3} = G_53; G{5,4} = G_54; G{5,5} = G_55; G{5,6} = G_56;
G{6,1} = G_61; G{6,2} = G_62; G{6,3} = G_63; G{6,4} = G_64; G{6,5} = G_65; G{6,6} = G_66;
end


