function result = SimulateBurst(params,s_idx)
% This function simulates the dynamic evolution of E-P distance and transcriptional burst.

%% Initialization
%% pre-allocating transcriptional variable and Initialization
transcription_state = zeros(1, params.simulation_reaction_step); % recording the state after each reaction
transcription_time  = zeros(1, params.simulation_reaction_step); % recording the time at each reaction
pol_at_promoter     = zeros(1,params.simulation_reaction_step); % recording when pol is at promoter
pol_at_elongation   = zeros(1,params.simulation_reaction_step); % recording when pol starts elongating
burst_duration_index = [0;0]; %
mRNA = 0;
tau_array   = zeros(1, params.simulation_reaction_step); 
transcription_state(1,1) = 1;
transcription_time(1,1)  = 0;
%% pre-allocating enhancer-promoter variable and Initialization
if params.EP_flag
    current_position = zeros(params.beads_num, 3);
    rand_array       = params.b.*randn(params.beads_num,params.dimension);
    for ps_idx = 2:params.beads_num
        current_position(ps_idx,:) = current_position(ps_idx-1,:) + rand_array(ps_idx,:);
    end
    % initialize rouse chain, harmonic potential between consecutive monomers
    interaction_matrix = InitializeConnectivityMatrix(params);
    % an equilibratio of 1e4 steps was run before taking samples,
    randn_array = randn(params.beads_num,params.dimension,1000000).*params.factor;
    for index = 1:1000000
        current_position = current_position -(1./params.friction_coef).*...
            interaction_matrix*(current_position*params.dt) + randn_array(:,:,index); 
    end
    EP_distance = pdist([current_position(params.enhancer_index,:);...
        current_position(params.promoter_index,:)]);
    i = 1;
    randn_array = randn(params.beads_num,params.dimension,100000).*params.factor;
end
%% Initialization end
%% Main simulation loop
tic;
current_idx = 2; record_steps = 0; 
while current_idx <= params.simulation_reaction_step && transcription_time(1,current_idx - 1) < params.simulation_time
    if params.EP_flag
        u1 = rand(1); H = 1; steps = 0; 
        while H > u1
            kon1 =  (EP_distance <= params.distance_T)*params.k_on1_max + ...
                (EP_distance > params.distance_T)*(params.k_on1 + (params.k_on1_max...
                - params.k_on1)/(1+((EP_distance-params.distance_T)/(params.distance_05...
                -params.distance_T))^params.H));
            krec =  (EP_distance <= params.distance_T)*params.k_recruitment_max + ...
                (EP_distance > params.distance_T)*(params.k_recruitment + (params.k_recruitment_max...
                - params.k_recruitment)/(1+((EP_distance-params.distance_T)/(params.distance_05...
                -params.distance_T))^params.H));
            krel =  (EP_distance <= params.distance_T)*params.k_release_max + ...
                (EP_distance > params.distance_T)*(params.k_release + (params.k_release_max...
                - params.k_release)/(1+((EP_distance-params.distance_T)/(params.distance_05...
                -params.distance_T))^params.H));
            kon2 = params.k_on2; koff2 = params.k_off2;koff1 = params.k_off1;koff3 = params.k_off3;
            if transcription_state(1,current_idx - 1) == 1; a_tot = kon2;  cutoffs = [0 kon2 kon2 kon2 kon2 kon2 kon2 kon2]; end
            if transcription_state(1,current_idx - 1) == 2; a_tot = kon1 + koff2; cutoffs = [0 0 kon1 kon1 + koff2 kon1 + koff2 kon1 + koff2 kon1 + koff2 kon1 + koff2]; end
            if transcription_state(1,current_idx - 1) == 3; a_tot = koff1 + krec; cutoffs = [0 0 0 0 koff1  koff1 + krec  koff1 + krec  koff1 + krec ]; end
            if transcription_state(1,current_idx - 1) == 4; a_tot = krel + koff3; cutoffs = [0 0 0 0 0 0 krel krel + koff3]; end
            H = H - a_tot*H*params.dt;
            current_position = current_position - (1./params.friction_coef).*...
                interaction_matrix*(current_position*params.dt) + randn_array(:,:,i);
            steps = steps + 1; i = i + 1; record_steps = record_steps + 1;
            if i == 100000; i = 1;randn_array = randn(params.beads_num,params.dimension,100000).*params.factor; end
            EP_distance = pdist([current_position(params.enhancer_index,:);current_position(params.promoter_index,:)]);            
        end
        tau = steps*params.dt;
    else
        kon1 = params.k_on1; krec = params.k_recruitment; krel = params.k_release;
        kon2 = params.k_on2; koff2 = params.k_off2;koff1 = params.k_off1;koff3 = params.k_off3;
        if transcription_state(1,current_idx - 1) == 1; a_tot = kon2;  cutoffs = [0 kon2 kon2 kon2 kon2 kon2 kon2 kon2]; end
        if transcription_state(1,current_idx - 1) == 2; a_tot = kon1 + koff2; cutoffs = [0 0 kon1 kon1 + koff2 kon1 + koff2 kon1 + koff2 kon1 + koff2 kon1 + koff2]; end
        if transcription_state(1,current_idx - 1) == 3; a_tot = koff1 + krec; cutoffs = [0 0 0 0 koff1 koff1 + krec  koff1 + krec  koff1 + krec ]; end
        if transcription_state(1,current_idx - 1) == 4; a_tot = krel + koff3; cutoffs = [0 0 0 0 0 0 krel krel + koff3]; end 
        u1 = rand(1);
        tau = roundn((1/a_tot)*log(1/u1),-2);
    end
    u2 = rand(1)*a_tot;
    reaction_list = [0 2 3 1 2 4 3 2];
    mu = reaction_list(find(u2<cutoffs,1,'first'));
    transcription_time(1,current_idx) = transcription_time(1,current_idx - 1) + tau;
    transcription_state(1,current_idx) = mu;
    tau_array(1,current_idx) = tau;
    if mu == 3 && transcription_state(1,current_idx - 1) == 4; pol_at_elongation(1,current_idx) = 1; end % if it does state 3-->2, that is elongation (1 mRNA produced)
    if mu == 4; pol_at_promoter(1,current_idx) = 1; end % if it does state 3, that is pol II at promoter
    if mu == 3 && transcription_state(1,current_idx - 1) == 2; burst_duration_index(1,end + 1) = current_idx; end
    if ((mu == 2) && (transcription_state(1,current_idx - 1) == 3 || transcription_state(1,current_idx - 1) == 4 )); burst_duration_index(2,end) = current_idx; end
    current_idx = current_idx + 1;
end % end loops
 timerVal = toc;
 X = ['Time taken for simulation',num2str(s_idx),': ',num2str(timerVal)];
 disp(X)
%% burst duration & mRNA
burst_duration_index = burst_duration_index(:,2:end);
if burst_duration_index(2,end) == 0; burst_duration_index = burst_duration_index(:,1:end - 1); end
for idx = 1:size(burst_duration_index,2)
    mRNA(end+1) = sum(pol_at_elongation(burst_duration_index(1,idx):burst_duration_index(2,idx)));
end
mRNA = mRNA(2:end);
%%
%% Save simulation data
record_size = params.simulation_reaction_step;
if current_idx < params.simulation_reaction_step; record_size = current_idx - 1; end
result.params              = params;
result.record_size         = record_size;
result.transcription_state = transcription_state(1:record_size);
result.transcription_time  = transcription_time(1:record_size);
result.pol_at_promoter     = pol_at_promoter(1:record_size);
result.pol_at_elongation   = pol_at_elongation(1:record_size);
result.mRNA                = mRNA;
result.tau_array           = tau_array(1:record_size);
result.burst_duration_index = burst_duration_index;
filename = sprintf('//%d.mat', s_idx);
save([params.result_base_folder,params.filename,filename],'result');
end





