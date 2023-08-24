

% main offline

tic
% the preprocessed input of concentration tensor has been calculated 
% or calculated by FFTbasedMicromech3D
clear all;
current_loc = pwd;
% Assign cluster number 
% Nm = 64;
% Ni = 8;
Nm = 64;
Ni1 = 32;
Ni2 = 32;
Ni3 = 32;
Ni4 = 32;
Ntotal = [Nm, Ni1, Ni2, Ni3, Ni4];
input_filename = 'multiphase_input.txt';
data_origin = 2; 
% 1 represents the inner data calculated by FFTbasedMicromech3D
% 2 represents the outer dataset
preprocess_micro_SCA(current_loc, Ntotal, input_filename, data_origin);

fprintf('finish offline stage. \n');
fprintf('cluster num: matrix %d, inclusion %d.',Ntotal(1), Ntotal(2:end));
%fprintf('cluster num: matrix %d, inclusion %d.\n',Nm, Ni1,Ni2,Ni3,Ni4);
toc