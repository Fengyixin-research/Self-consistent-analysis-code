function [C] = C0_cal_3D(lamda,mu)
C = zeros(6,6);
C(1:3,1:3) = lamda;
C(1,1) = lamda + 2*mu;
C(2,2) = C(1,1); C(3,3) = C(1,1);
C(4,4) = mu; 
C(5,5) = C(4,4); C(6,6) = C(4,4);
end