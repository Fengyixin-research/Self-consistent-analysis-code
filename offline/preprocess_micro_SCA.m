function preprocess_micro_SCA(current_loc, Ntotal,input_filename, data_origin)
% number of material phases
N_phase = length(Ntotal);
Nt = sum(Ntotal); 
Nm = Ntotal(1); % first is assumed with matrix
Ni_total = Ntotal(2:end); % summerize the inclusion phases
%Nt = Nm + Ni1 + Ni2 + Ni3 + Ni4;
% preprocess of data
[phase,material,C0,location,prestrain] = preprocess(input_filename, data_origin);

% K_means clustering
[matrixlist,inclusionlist_total,mech_int_part1, mech_int_part2,vf,idx] = patternA_3D(N_phase, Nt, Nm, Ni_total, phase, material, C0, data_origin, prestrain);

% store stffness tensor
C_total = cell(N_phase, 1);
for i = 1:N_phase
    C_total{i,1} = C_cal_3D_orthotropic_comsol(material(i,1), material(i,2), material(i,3), material(i,4), material(i,5), material(i,6), material(i,7), material(i,8), material(i,9));
end

save('preprocess.mat','Nt','Nm','Ni_total','C0','mech_int_part1',...
    'mech_int_part2','vf','C_total','material');

save('post_location.mat', 'location','idx','matrixlist','inclusionlist_total','Ni_total','Nm');

end

%% preprocess 
function [phase,material,C0,location,prestrain] = preprocess(input_filename, data_origin) %main
% orthotropic,E1,E2,E3,v12,v23,v13,G12,G23,G13
% m,n,p:number of points along direction x,y,z
% material:lamda,mu
% phase 
% C0

% filename = 'stress_strain_2dplastic.txt';
data = load(input_filename,'-regexp','^(?/%)...');
x=data(:,1); y=data(:,2); z=data(:,3); material_property=data(:,4); 
prestrain = data(:, 5:40); 
E1 = data(:,41); E2 = data(:,42); E3 = data(:,43);
v12 = data(:,44); v23 = data(:,45); v13 = data(:,46);
G12 = data(:,47); G23 = data(:,48); G13 = data(:,49);

% location of points
location = [x, y, z];

x_u = unique(x);
y_u = unique(y);
z_u = unique(z);
m = round((max(x_u) - min(x_u))/(x_u(2)-x_u(1))) + 1;
n = round((max(y_u) - min(y_u))/(y_u(2)-y_u(1))) + 1;
p = round((max(z_u) - min(z_u))/(z_u(2)-z_u(1))) + 1;
phase = reshape(material_property, m, n, p); 

material_property_unique = unique(material_property);
E_total = zeros(size(material_property_unique, 1), 3);
nu_total = zeros(size(material_property_unique, 1), 3);
mu_total = zeros(size(material_property_unique, 1), 3);
for i = 1:length(material_property_unique)
    E_total(i, :) = [mean(E1(material_property==material_property_unique(i,1))), mean(E2(material_property==material_property_unique(i,1))), mean(E3(material_property==material_property_unique(i,1)))];
    nu_total(i, :) = [mean(v12(material_property==material_property_unique(i,1))), mean(v23(material_property==material_property_unique(i,1))), mean(v13(material_property==material_property_unique(i,1)))];
    mu_total(i, :) = [mean(G12(material_property==material_property_unique(i,1))), mean(G23(material_property==material_property_unique(i,1))), mean(G13(material_property==material_property_unique(i,1)))];
    if i == length(material_property_unique)  % the prescribed E and nu for copper, which is not a orthotropic material
        E_total(i, :) = 1.2758e11;
        nu_total(i, :) = 0.347;
        mu_total(i, :) = 1.2758e11/(2*(1+0.347));
    end
end
material = [E_total,nu_total,mu_total,material_property_unique];

E0 = 1/2*(max(material(:,1:3),[],'all')+min(material(:,1:3),[],'all'));
nu0 = 1/2*(max(material(:,4:6),[],'all')+min(material(:,4:6),[],'all'));

lamda0 = E0*nu0/((1+nu0)*(1-2*nu0));
mu0 = E0/(2*(1+nu0));
C0 = C0_cal_3D(lamda0,mu0);
end

%% compute interaction and volume_fraction
function [matrixlist,inclusionlist_total,mech_int_part1, mech_int_part2, vf,idx] = patternA_3D(N_phase, Nt, Nm, Ni_total, phase, material, C0, data_origin, prestrain)
% Input:
% Nm, Ni, phase: 0/1, material: [lamda,mu, material_property], C0
%

% size of the phase
[m,n,p] = size(phase);

%% Generate the gvector
% Discrete frequency point, assuming the length of each element is 1
em = (ceil(-m/2):floor((m-1)/2));
en = (ceil(-n/2):floor((n-1)/2));
ep = (ceil(-p/2):floor((p-1)/2));
display('Start assembling gvector ...');

% Separate the elements in the inclusion and matrix domains
meshvector = reshape(phase,[],1);
matrixlist = meshvector == N_phase-1;
inclusionlist_total = cell(N_phase-1, 1);
for i = 1:N_phase-1
    inclusionlist_total{i, 1} = meshvector == i-1;
end

data_loadtime = 2; 

tic
if data_origin ~= 2 % 
    [e11,e22,e33,e12,e23,e13] = FFTbasedMicromech3D(phase, material, C0, [1,0,0,0,0,0]');
    A = [reshape(e11,[],1),reshape(e22,[],1),reshape(e33,[],1),reshape(e12,[],1),reshape(e23,[],1),reshape(e13,[],1)];
    [e11,e22,e33,e12,e23,e13] = FFTbasedMicromech3D(phase, material, C0, [0,1,0,0,0,0]');
    A = [A,reshape(e11,[],1),reshape(e22,[],1),reshape(e33,[],1),reshape(e12,[],1),reshape(e23,[],1),reshape(e13,[],1)];
    [e11,e22,e33,e12,e23,e13] = FFTbasedMicromech3D(phase, material, C0, [0,0,1,0,0,0]');
    A = [A,reshape(e11,[],1),reshape(e22,[],1),reshape(e33,[],1),reshape(e12,[],1),reshape(e23,[],1),reshape(e13,[],1)];
    [e11,e22,e33,e12,e23,e13] = FFTbasedMicromech3D(phase, material, C0, [0,0,0,1,0,0]');
    A = [A,reshape(e11,[],1),reshape(e22,[],1),reshape(e33,[],1),reshape(e12,[],1),reshape(e23,[],1),reshape(e13,[],1)];
    [e11,e22,e33,e12,e23,e13] = FFTbasedMicromech3D(phase, material, C0, [0,0,0,0,1,0]');
    A = [A,reshape(e11,[],1),reshape(e22,[],1),reshape(e33,[],1),reshape(e12,[],1),reshape(e23,[],1),reshape(e13,[],1)];
    [e11,e22,e33,e12,e23,e13] = FFTbasedMicromech3D(phase, material, C0, [0,0,0,0,0,1]');
    A = [A,reshape(e11,[],1),reshape(e22,[],1),reshape(e33,[],1),reshape(e12,[],1),reshape(e23,[],1),reshape(e13,[],1)];
else
    if data_loadtime == 2
        A = prestrain;
    end
end
toc

% Divide the gvector according to the material
gvector_m = A(matrixlist,:);
gvector_i_total = cell(N_phase-1, 1);
for i = 1:N_phase-1
    gvector_i_total{i,1} = A(inclusionlist_total{i,1},:); 
end

% Save the gvector data to the file
display('Finish assembling gvector');

%% k-mean clustering based on the conectration tensor (gvector)
display('Start kmeans clustering ...');

idx_m = kmeans(gvector_m,Nm,'MaxIter',5000);
idx_i_total = cell(N_phase-1, 1);
for i = 1:N_phase-1
    idx_i_total{i,1} = kmeans(gvector_i_total{i,1},Ni_total(i),'MaxIter',5000);
end

% The total index of the patterns
idx = zeros(m*n*p,1);
idx(matrixlist) = idx_m;
for i = 1:N_phase-1
    idx(inclusionlist_total{i,1}) = idx_i_total{i,1};
end

% compute the interaction matrice
display('Start computing the interaction matrice...');
[mech_int_part1, mech_int_part2] = SCA_mech_int_cal(m, n, p, em, en, ep, inclusionlist_total, matrixlist, Ni_total, Nm, idx_i_total, idx_m, idx);

%% compute the volume fraction
%Nt = Nm + sum(Ni_total);
vf = zeros(Nt, 1);
for i = 1:N_phase
    if i ~= N_phase
        current_idx = idx_i_total{i,1};
        current_N = Ni_total(i);
        min_N = sum(Ni_total(1:i-1));
    else
        current_idx = idx_m;
        current_N = Nm;
        min_N = sum(Ni_total);
    end
    for j = 1:current_N
        vf(min_N+j) = sum(current_idx==(j))/m/n/p;
    end
end
end
		
