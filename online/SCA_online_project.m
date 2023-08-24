clear all;

current_loc = pwd;
load('preprocess.mat');

full_strain = [0,0,0.007,0,0,0];
count_total = 10;
disp_factor = linspace(0,1,count_total+1);
disp_fac = zeros(1, length(disp_factor)-1);
for i = 1:length(disp_factor)-1
    disp_fac(i) = disp_factor(i+1) - disp_factor(i);
end

H_table = [0.0372e9 1.1887e9 0.7943];

%% Initialization
etotal = zeros(6*Nt,1);
sigma = zeros(6*Nt,1);
eqpl = zeros(Nt,1);

E_initial = sum(reshape(etotal(1:6*Nt),6,Nt)*diag(vf),2);
Sigma_initial = sum(reshape(sigma,6,Nt)*diag(vf),2);

NewtonGol = 100; 
FixE0 = 50; 
toler = 1e-3;
NewtonUmat = 100;

etotal_initial = etotal;
sigma_initial = sigma;
eqpl_initial = eqpl;
C0_initial = C0;
lamda0_initial = C0_initial(2,1);
mu0_initial = C0_initial(4,4);
E0_initial = (lamda0_initial*3+mu0_initial*2)*mu0_initial/(lamda0_initial+mu0_initial);
v0_initial = lamda0_initial/(2*(lamda0_initial+mu0_initial));

Stress_out = zeros(6, count_total);
Strain_out = zeros(6, count_total);

J = zeros(6,6);
J(1:3,1:3) = 1/3;
K = zeros(6,6);
K(1:3,1:3) = -1/3;
K(1,1) = 2/3; K(2,2) = 2/3; K(3,3) = 2/3;
K(4,4) = 1/2; K(5,5) = 1/2; K(6,6) = 1/2;

% loading
for total_iiter = 1:count_total
    display(['Start calculate ',num2str(total_iiter),' iiter.']);
    dE = zeros(6,1);
    dE = disp_fac(total_iiter)*full_strain' + dE;

    etotal = etotal_initial;
    sigma = sigma_initial;
    eqpl = eqpl_initial;
    C0 = C0_initial;
    lamda0 = lamda0_initial;
    detotal = zeros(6*(Nt),1);    

    error = 0;
    iterE0 = 0;
    deltalamda0 = lamda0;

    %% while loop  
    while abs(deltalamda0/lamda0) > toler && iterE0 < FixE0
        iterE0 = iterE0 + 1;
        lamda0 = C0(1,2); mu0 = C0(4,4);
        c1 = 1/4/mu0; c2 = (lamda0+mu0)/mu0/(lamda0+2*mu0);  
        S = c1*mech_int_part1+c2*mech_int_part2;  

        Cdiag = cell(1, Nt);
        for i = 1:Nt
            [Cdiag{i}] = C0;
        end
        S_Ct = S*sparse(blkdiag(Cdiag{:}));        
        dsigma = zeros(6*Nt,1);    
        iteration = 0;
        % Newton iteration
        for i = 1:NewtonGol
            dsigma = reshape(dsigma,6,Nt);
            dCalgdiag = cell(1,Nt);
            
            for iiter_phase = 1:length(Ni_total)+1
                if iiter_phase <= length(Ni_total)
                    current_N = Ni_total(iiter_phase);
                    current_C = C_total{iiter_phase, 1};
                    for j = 1:current_N
                        current_j = j + sum(Ni_total(1:iiter_phase-1));
                        [Ddsdde_j,~,~,dsigma_i] = umatelastic_3D(detotal(6*current_j-5:6*current_j),...
                            eqpl_initial(current_j),sigma_initial(6*current_j-5:6*current_j),current_C);     
                        [dCalgdiag{current_j}] = Ddsdde_j-C0;
                        dsigma(:,current_j) = dsigma_i;
                    end
                else
                    current_N = Nm;
                    current_C = C_total{iiter_phase, 1};
                    for j = 1:current_N
                        current_j = j + sum(Ni_total);
                        [Ddsdde_j,~,~,dsigma_i] = umatJ2_3D_analy(detotal(6*current_j-5:6*current_j),...
                            eqpl_initial(current_j),sigma_initial(6*current_j-5:6*current_j), current_C(2,1), current_C(4,4), H_table, NewtonUmat);    
                        [dCalgdiag{current_j}] = Ddsdde_j-C0;            
                        dsigma(:,current_j) = dsigma_i;
                    end
                end
            end
            dsigma = reshape(dsigma,[],1);
            dCalg = sparse(blkdiag(dCalgdiag{:}));

            % residual
            r(1:6*Nt,1) = detotal(1:6*Nt,1) + S*dsigma - S_Ct*(detotal(1:6*Nt,1)) - repmat(dE,Nt,1);        
            % Jacobian
            M = zeros(6*(Nt),6*(Nt));
            M(1:6*Nt,1:6*Nt) = S*dCalg + eye(6*Nt);
            % Compute the new detotal
            Ddetotal = -M\r;
            detotal = detotal + Ddetotal;         

            iteration = iteration + 1;
            re = zeros(Nt,1);
            for m = 1:Nt
                re(m) = norm(Ddetotal((6*m-5):6*m))/norm(detotal((6*m-5):6*m));
            end        
            if max(re)< 1e-8
                break;
            end        
        end
        % Compute the effective modulus of the material
        Minv = inv(M);
        Ceff = zeros(6,6);
        for j = 1:Nt
            Ceff = Ceff + vf(j)*(dCalgdiag{j}+C0)*(Minv(6*j-5:6*j,:)*repmat(eye(6),Nt,1));
        end
        % projection
        C0 = sum(J.*Ceff,'all')*J + 1/5*sum(K.*Ceff,'all')*K;    
        deltalamda0 = (C0(1,2) - lamda0);
        lamda0 = C0(1,2);
        if iterE0 == FixE0
            display(['No convergence of E0: ',num2str(abs(deltalamda0/lamda0))]);
        end    
    end


    % update 
    for iiter_phase = 1:length(Ni_total)+1
        if iiter_phase <= length(Ni_total)
            current_N = Ni_total(iiter_phase);
            current_C = C_total{iiter_phase, 1};
            for j = 1:current_N
                current_j = j + sum(Ni_total(1:iiter_phase-1));
                [~,eqplj,sigmaj,~] = umatelastic_3D(detotal(6*current_j-5:6*current_j),...
                    eqpl_initial(current_j),sigma_initial(6*current_j-5:6*current_j),current_C);             
                etotal(6*current_j-5:6*current_j) = etotal(6*current_j-5:6*current_j) + detotal(6*current_j-5:6*current_j);
                eqpl(current_j) = eqplj;
                sigma(6*current_j-5:6*current_j) = sigmaj;     
            end
        else
            current_N = Nm;
            current_C = C_total{iiter_phase, 1};
            for j = 1:current_N
                current_j = j + sum(Ni_total);
                [~,eqplj,sigmaj,~] = umatJ2_3D_analy(detotal(6*current_j-5:6*current_j),...
                eqpl_initial(current_j),sigma_initial(6*current_j-5:6*current_j), current_C(2,1), current_C(4,4), H_table, NewtonUmat);        
                etotal(6*current_j-5:6*current_j) = etotal(6*current_j-5:6*current_j) + detotal(6*current_j-5:6*current_j);
                eqpl(current_j) = eqplj;
                sigma(6*current_j-5:6*current_j) = sigmaj;   
            end
        end
    end
    
    E = sum(reshape(etotal,6,Nt)*diag(vf),2);
    Sigma = sum(reshape(sigma,6,Nt)*diag(vf),2);
    etotal_initial = etotal;
    sigma_initial = sigma;
    eqpl_initial = eqpl;
    C0_initial = C0;
    lamda0_initial = lamda0;   
    Strain_out(:, total_iiter) = E;
    Stress_out(:, total_iiter) = Sigma;     
end