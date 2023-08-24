function [ddsdde_out, eqpl_out, sigma_out, dsigma_out] =... 
umatJ2_3D_analy(detotal, eqpl, sigma, lamda, mu, H_table, NewtonUmat)
% elastoplastic material law
% input
% detotal: total strain increment
% eqpl: existed equivalent plastic strain
% sigma: existed sigma
% E: Young's modulus
% v: Poisson's ratio
% H_table: f(x) = a + b*x.^c, a:initial yield stress, b:hard coefficient,
% c:index, x: equivalent plastic strain
% NewtonUmat: maxiteration

% output
% ddsdde_out: algorithm tensor
% eqpl_out: equivalent plastic strain 
% sigma_out
% dsigma_out: stress increment

ebulk3 = 3*lamda + 2*mu; % bulk modulus*3
% elastic
ddsdde = zeros(6,6);
for i = 1:6
    if i < 4
        for j = 1:3
            ddsdde(i,j) = lamda;
        end
        ddsdde(i,i) = lamda + 2*mu;
    else
        ddsdde(i,i) = mu;
    end
end

% attempt stress
dsigmatotal = ddsdde*detotal; 
sigma0 = sigma;
sigma = sigma + dsigmatotal; 

% Mises stress
smises = sqrt(((sigma(1)-sigma(2))^2 + (sigma(2)-sigma(3))^2 + (sigma(3)-sigma(1))^2 + ...
    6*norm(sigma(4:6))^2)/2);
% hard, syield
[hard, syield] = hardsub(eqpl, H_table);
%syield0 = syield;
% Mises stress> yield stress
toler = 1e-7;
if (smises > (1+toler)*syield) 
    shydro = sum(sigma(1:3))/3; 
    flow = zeros(6,1);
    flow(1:3) = (sigma(1:3) - shydro)/smises;
    flow(4:6) = sigma(4:6)/smises;
    
    deqpl = 0;
    for i = 1:NewtonUmat
        rhs = smises - 3*mu*deqpl - syield;
        deqpl = deqpl + rhs/(3*mu + hard);
        [hard, syield] = hardsub(eqpl+deqpl, H_table);
        if (abs(rhs) < toler*syield); break; end
    end
     
    sigma_out = zeros(6,1);
    for i = 1:6
        if i < 4
            sigma_out(i) = flow(i)*syield + shydro;
        else
            sigma_out(i) = flow(i)*syield;  
        end
    end
    eqpl_out = eqpl + deqpl;
    dsigma_out = sigma_out - sigma0; 
    
    % Jacobian
    effmu = mu*syield/smises;
    efflamda = (ebulk3 - 2*effmu)/3;
    effhard = 3*mu*hard/(3*mu+hard) - 3*effmu;
    ddsdde_out = zeros(6,6);
    for i = 1:6
        if i < 4
            for j = 1:3
                ddsdde_out(i,j) = efflamda;
            end
            ddsdde_out(i,i) = efflamda + 2*effmu;
        else
            ddsdde_out(i,i) = effmu;
        end
    end
    for i = 1:6
        for j = 1:6
            ddsdde_out(i,j) = ddsdde_out(i,j) + effhard*flow(i)*flow(j);
        end
    end
    
else
    ddsdde_out = ddsdde;
    eqpl_out = eqpl;
    sigma_out = sigma;
    dsigma_out = dsigmatotal;
end
end

function [hard_out, syield_out] = hardsub(eqpl,H_table)
syield0 = H_table(1);
H = H_table(2);
n = H_table(3);

hard_out = H*(n*(eqpl+1e-5)^(n-1));
syield_out = syield0 + H*((eqpl+1e-5)^n -1e-5^n);

end





