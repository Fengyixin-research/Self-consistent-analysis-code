function [D] = C_cal_3D_orthotropic_comsol(E1,E2,E3,nu12,nu23,nu13,mu12,mu23,mu13)
% stiffness tensor for orthotropic
% E1,E2,E3,nu12,nu23,nu13,mu12,mu23,mu13

D = zeros(6,6);
D_denom = E2*E3*nu13^2 - E1*E2 + 2*nu12*nu23*nu13*E2*E3 + E1*E3*nu23^2 + E2^2*nu12^2;
D(1,1) = E1^2*(E3*nu23^2-E2)/D_denom; D(1,2) = -E1*E2*(E3*nu23*nu13+E2*nu12)/D_denom;
D(1,3) = -E1*E2*E3*(nu12*nu23 + nu13)/D_denom; D(2,2) = E2^2*(E3*nu13^2 - E1)/D_denom;
D(2,3) = -E2*E3*(E2*nu12*nu13 + E1*nu23)/D_denom; D(3,3) = E2*E3*(E2*nu12^2 - E1)/D_denom;

D(2,1) = D(1,2); D(3,1) = D(1,3); D(3,2) = D(2,3);

D(4,4) = mu12; D(5,5) = mu23; D(6,6) = mu13;

end