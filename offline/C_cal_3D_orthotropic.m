function [C] = C_cal_3D_orthotropic(E1,E2,E3,nu12,nu23,nu13,mu12,mu23,mu13)
% stiffness tensor for orthotropic
% E1,E2,E3,nu12,nu23,nu13,mu12,mu23,mu13

S = zeros(6,6);
S(1:3,1:3) = [1/E1        -nu12/E2    -nu13/E3;
              -nu12/E2    1/E2        -nu23/E3;
              -nu13/E1    -nu23/E2    1/E3;];
S(4,4) = 1/mu12;
S(5,5) = 1/mu23;
S(6,6) = 1/mu13;

C = inv(S);
end