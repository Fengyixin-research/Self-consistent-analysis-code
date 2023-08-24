function [ddsdde_out, eqpl_out, sigma_out, dsigma_out] =... 
umatelastic_3D(detotal, eqpl, sigma, C)
% linear elastic
ddsdde = C;

% 
dsigmatotal = ddsdde*detotal; 
sigma = sigma + dsigmatotal; 

ddsdde_out = ddsdde;
eqpl_out = eqpl;
sigma_out = sigma;
dsigma_out = dsigmatotal;
end