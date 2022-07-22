clear;clc;
close all;

K = [0.1,0.2];
for idx = 1:size(K,2)
    % for idx = 1:size(enhancer,2)
    %% Parameters setting
    % Upstream: Chromatin conformation
    input_options.EP_flag           = true;
    input_options.attraction_coef   = K(idx); % K, E-P communication coefficient
    input_options.enhancer_index    = 25; % enhancer index number
    input_options.promoter_index    = 75; % promoter index number
    
    input_options.simulation_time   = 100000; % [s]
    input_options.result_base_folder = fullfile(pwd, 'Results');
    
    %% NumericalSimulation
    EPDistNumerical(input_options);

    %% Theroetical
    [EP_Bin_Thero,EPPDF_Thero] = EPDistTheroetical(input_options);
    EPPDF.K(idx).EP_Bin_Thero = EP_Bin_Thero;
    EPPDF.K(idx).EPPDF_Thero = EPPDF_Thero;
end
save EPPDF

%% draw figures Figure S10A
figure1 = figure;
load("Results\E_25_P_75_KEP_1.000000e-01.mat")
h = histogram(distance(1,:),30,'Normalization','pdf');
hold on
plot(EPPDF.K(1).EP_Bin_Thero,EPPDF.K(1).EPPDF_Thero);
set(figure1,'position',[300 400 280 190]);
axis([-0.01 1 0 6])

figure2 = figure;
load("Results\E_25_P_75_KEP_2.000000e-01.mat")
h = histogram(distance(1,:),25,'Normalization','pdf');
hold on
plot(EPPDF.K(2).EP_Bin_Thero,EPPDF.K(2).EPPDF_Thero);
set(figure2,'position',[300 400 280 190]);
axis([-0.01 1 0 6])



function EPDistNumerical(input_options)
params = ParametersBurst(input_options);
distance  = zeros(1,params.simulation_time);
current_position = zeros(params.beads_num, 3);
rand_array       = params.b.*randn(params.beads_num,params.dimension);
for ps_idx = 2:params.beads_num
    current_position(ps_idx,:) = current_position(ps_idx-1,:) + rand_array(ps_idx,:);
end
% Initialization
Rouse_matrix = InitializeConnectivityMatrix(params);
for i = 1:10000
    current_position = current_position - (1./params.friction_coef).*...
        Rouse_matrix*(current_position*params.dt) + randn(params.beads_num,params.dimension).*params.factor;
end

d_idx = 1;
for i = 1:params.simulation_time/params.dt
    current_position = current_position - (1./params.friction_coef).*...
        Rouse_matrix*(current_position*params.dt)+ randn(params.beads_num,params.dimension).*params.factor;
    if mod(i,100) == 0
        current_distance = [current_position(params.enhancer_index,:);current_position(params.promoter_index,:)];
        distance(d_idx) = pdist(current_distance);
        d_idx = d_idx + 1;
    end
end
filename = sprintf('//E_%d_P_%d_KEP_%d.mat',params.enhancer_index,params.promoter_index,params.attraction_coef);
save([params.result_base_folder,filename],'distance');
end


function [EP_Bin_Thero,EPPDF_Thero] = EPDistTheroetical(input_options)
params = ParametersBurst(input_options);
Kaa = sqrt(params.diffusion_const*params.friction_coef*...
    (params.spring_const(1,1)/(params.promoter_index - params.enhancer_index)...
    +params.attraction_coef)^(-1));
dt = 0.01;
d_EP = 0.01:dt:5;
EP_Bin_Thero = d_EP;
EPPDF_Thero = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
end




