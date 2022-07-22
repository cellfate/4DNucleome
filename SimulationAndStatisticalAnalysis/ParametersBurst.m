function params = ParametersBurst(input_options)
params.simulated_on   = true;
params.EP_flag        = true;
params.compute_TR     = true;
params.dimension      = 3; % spatial dimension
params.beads_num      = 100; % beads number
params.enhancer_index = 25; % enhancer index number
params.promoter_index = 75; % promoter index number
params.dt             = 0.01;  % time step
params.friction_coef  = 10; % friction coefficient [Newton][sec]/[mu m]
params.diffusion_const = 4e-3/params.friction_coef; % diffusion constant [mu m]^2 /[sec]
params.b              = 0.1; % STD of distance between adjacent monomers
params.spring_const   = ones(params.beads_num); % spring constant
params.attraction_coef = 1; % attraction interaction constant
params.factor         = sqrt(2*params.diffusion_const*params.dt);
params.distance_T     = 0.1; 
params.distance_05    = 0.2;
params.H              = 3;
params.encounter_dist = params.distance_T; % interaction distance

params.k_on1_max         = 0.020; % burst initiation rate
params.k_on1             = 0.003; % burst initiation rate
params.k_recruitment_max = 0.5; % burst polymerase recruitment rate
params.k_recruitment     = 0.010; % burst polymerase recruitment rate
params.k_release_max     = 0.3; % burst polymerase pause release rate
params.k_release         = 0.008; % burst polymerase pause release rate
params.k_off2            = 0.008;% off1 to off2
params.k_on2             = 0.008;% off2 to off1
params.k_off1            = 0.003;% rec to off1
params.k_off3            = 0.003;% rec to off2

params.simulation_num    = 100; % simulation number
params.simulation_reaction_step = 5000;
params.simulation_time   = 10000;
params.snapshot_interval = 20; % [s]

params.result_base_folder = fullfile(pwd,'Results');
params.filename           = [];
params.elongation_time    = 100;

%% Read the input options
intput_params = fieldnames(input_options);
for i = 1:length(intput_params)
    name  = intput_params{i};
    value = input_options.(name);
    if isfield(params, name)
        params.(name) = value;
    end
end


