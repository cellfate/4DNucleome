function result = AnalyseBurst(params)
if params.simulated_on
    % This function analyzes data (numerical & theroetical) after parallel computing
    %% =====================numerical =======================================
    transcription_state = zeros(params.simulation_num,params.simulation_reaction_step);
    transcription_time = zeros(params.simulation_num,params.simulation_reaction_step);
    pol_at_elongation = zeros(params.simulation_num,params.simulation_reaction_step);
    record_size = params.simulation_reaction_step;
    snap_transcription_state = zeros(params.simulation_num, 10000);
    mRNA = [];
    snap_size = 10000;
    total_on_time = [];
    total_off_time = [];
    %% load data
    for i = 1:1:params.simulation_num
        filename = sprintf('//%d.mat',i);
        temp = load([params.result_base_folder,params.filename,filename]);
        burst_duration_index = temp.result.burst_duration_index;
        mRNA = [mRNA,temp.result.mRNA];
        if record_size > temp.result.record_size; record_size = temp.result.record_size; end
        transcription_state(i,1:record_size) = temp.result.transcription_state(1:record_size);
        transcription_time(i,1:record_size) = temp.result.transcription_time(1:record_size);
        pol_at_elongation(i,1:record_size) = temp.result.pol_at_elongation(1:record_size);
        result.simulate(i).data = temp;
        A = temp.result.transcription_time(burst_duration_index);
        total_on_time = [total_on_time,(A(2,1:end-1) - A(1,1:end-1))];
        total_off_time = [total_off_time,(A(1,2:end) - A(2,1:end-1))];

    end
    transcription_state = transcription_state(:,1:record_size);
    transcription_time =  transcription_time(:,1:record_size);
    pol_at_elongation = pol_at_elongation(:,1:record_size);
    if params.compute_TR
        %% Taking sample
        % each different simulation takes a different amount of time for each step transition (from 1 column to the next)
        % we now want to take equally spaced snapshots of the state of all simulations so we can visualize the data
        for i = 1:params.simulation_num
            temp_snap_size = floor(transcription_time(i,end)/params.snapshot_interval);
            % if the time for a designated snapshot(k) falls between col n and col n+1 of
            % one simulation (simulation m), record the n state as the state of the m simulation
            % at snapshot time k
            for k = 1:temp_snap_size
                index = find(transcription_time(i,:) > k*params.snapshot_interval, 1, 'first');
                snap_index = index - 1;
                snap_transcription_state(i, k) = transcription_state(i, snap_index);
            end
            if temp_snap_size < snap_size; snap_size = temp_snap_size; end
        end
        snap_transcription_state = snap_transcription_state(:,1:snap_size);
        
        %% calculating our predictions of experimental data
        % now we want to calculate the pol2 traveling ratio,
        % active transcribing fraction and nascent RNA Intensity
        %% each 'gene copy'(each simulation iteration) has two sections where polymerase can be
        %% 1) a promoter (which is bound by pol if a gene is in state 3)
        pol_at_pro_avg = sum(snap_transcription_state == 4)./params.simulation_num;
        % experiments in the literature have shown that only 1 polymerase can be
        % bound on 1 promoter at a time! so this makes sense
        
        %% 2) a gene body (which is bound by pol if a gene transitions from state 3-->2)
        % so that pol elongates and will fall off after a period of time ('reaches the end of the gene')
        elongation_time = params.elongation_time; %
        elong_start = pol_at_elongation.*transcription_time;
        elong_end = elong_start + (elong_start~=0)*elongation_time;
        % now we'll check which gene copies have a polymerase which is elongating in the gene body at snapshot k
        % elong_start_index = elong_start/params.snapshot_interval;
        % elong_end_index = elong_end/params.snapshot_interval;
        % elong_index = elong_end_index - elong_start_index;
        pol_at_elong = zeros(params.simulation_num,snap_size);
        for i = 1:params.simulation_num
            for j = 1:record_size
                for k = 1:snap_size
                    if k*params.snapshot_interval > elong_start(i,j) && k*params.snapshot_interval < elong_end(i,j)
                        pol_at_elong(i,k) = pol_at_elong(i,k) + 1;
                    end
                end
            end
        end
        pol_at_elong_avg = sum(pol_at_elong)/params.simulation_num;
        
        %% compute values for average pol2 occupancy and traveling ratio,
        % after simulations stabilize and before many of the m simulations start to end
        traveling_ratio_avg = pol_at_elong_avg./pol_at_pro_avg;
        pol_at_elong_avg(find(pol_at_elong_avg==0)) = [];
        pol_at_pro_avg(find(pol_at_pro_avg==0)) = [];
        traveling_ratio_avg(find(isnan(traveling_ratio_avg))) = [];
        traveling_ratio_avg(find(traveling_ratio_avg == inf)) = [];
    end
end
%% =====================theroteical =======================================
if params.EP_flag == false
    kon1 = params.k_on1; krec = params.k_recruitment; krel = params.k_release;
    kon2 = params.k_on2; koff2 = params.k_off2;koff1 = params.k_off1;koff3 = params.k_off3;
    % BS
    burst_size_theor = (krec.*krel./(krec.*koff3+(krel+koff3).*koff1));
    % ON
    total_on_time_theor = (krel+krec+koff3)./(krec*koff3+(krel+koff3)*koff1);
    % OFF
    total_off_time_theor = (koff2+kon2)./(kon2.*kon1);
    % CT
    total_cycle_time_theor = total_off_time_theor + total_on_time_theor;
    burst_frequency_theor = 1/total_cycle_time_theor ;
    traveling_ratio_theor = params.elongation_time*params.k_release;
else
    Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
        (abs(params.promoter_index - params.enhancer_index))+params.attraction_coef)^(-1));
    dt = 0.01;
    d_EP = 0.01:dt:5;
    P = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    
    d_T = params.distance_T;
    d_05 = params.distance_05;
    Hill = (d_EP > d_T).* (1./(1+((d_EP-d_T)./(d_05-d_T)).^params.H));
    Hill(1:d_T/dt) = 0;
    lambda_on1 = (d_EP <= d_T).*params.k_on1_max + (d_EP > d_T).*...
        (params.k_on1 + (params.k_on1_max -params.k_on1).*Hill);
    lambda_recruitment = (d_EP <= d_T).*params.k_recruitment_max + (d_EP > d_T).*...
        (params.k_recruitment + (params.k_recruitment_max -params.k_recruitment).*Hill);
    lambda_release = (d_EP <= d_T).*params.k_release_max + (d_EP > d_T).*...
        (params.k_release+ (params.k_release_max -params.k_release).*Hill);
    lambda_on2 = params.k_on2; lambda_off2 = params.k_off2;
    lambda_off1 = params.k_off1;lambda_off3 = params.k_off3;
    on1 = sum(lambda_on1.*P.*dt);
    recruitment = sum(lambda_recruitment.*P.*dt);
    release = sum(lambda_release.*P.*dt);
    on2 = sum(lambda_on2.*P.*dt);
    off2 = sum(lambda_off2.*P.*dt);
    off1 = sum(lambda_off1.*P.*dt);
    off3 = sum(lambda_off3.*P.*dt);
    % BS
    burst_size_theorV1 = (recruitment.*release./(recruitment.*off3+(release+off3).*off1));
    % ONPDF
    total_on_time_theorV1 = (recruitment+release+off3)/(recruitment.*off3+(release+off3).*off1);
    % OFF
    total_off_time_theorV1 = (off2+on2)./(on2.*on1);
    % tr
    traveling_ratio_theor = params.elongation_time*release;
    
    %% V2
    burst_size_theorV2 = sum((lambda_recruitment.*lambda_release./...
        (lambda_recruitment.*lambda_off3+...
        (lambda_release+lambda_off3).*lambda_off1)).*P.*dt);
    % ON
    total_on_time_theorV2 = sum((lambda_recruitment+lambda_release+lambda_off3)./...
        (lambda_recruitment.*lambda_off3+(lambda_release+lambda_off3).*lambda_off1).*P.*dt);
    % OFF
    total_off_time_theorV2 = sum((lambda_off2+lambda_on2)./(lambda_on2.*lambda_on1).*P.*dt);
    
    %% mu
    minlambda = min([params.k_on1,params.k_recruitment,params.k_release]);
    maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
        params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
        99,1, 'first'))/params.friction_coef;
    mu = minlambda/maxvelocity;
    burst_size_theor = 1./(1+mu).*burst_size_theorV1+mu./(1+mu).*burst_size_theorV2;
    total_on_time_theor = 1./(1+mu).*total_on_time_theorV1+mu./(1+mu).*total_on_time_theorV2;
    total_off_time_theor = 1./(1+mu).*total_off_time_theorV1+mu./(1+mu).*total_off_time_theorV2;
    burst_frequency_theor = 1./(total_on_time_theor+total_off_time_theor);
end


%% ====================Saving data=========================
if params.simulated_on
    result.snap_transcription_state = snap_transcription_state;
      if params.compute_TR; result.pol_at_elong = pol_at_elong; end
    result.mRNA = mRNA;
    result.total_on_time = total_on_time;
    result.total_off_time = total_off_time;
end


if params.simulated_on;result.burst_size = mean(mRNA);end
result.burst_size_theor = burst_size_theor;

if params.simulated_on;result.on_time = mean(total_on_time);end
result.on_time_theor = total_on_time_theor;

if params.simulated_on;result.off_time = mean(total_off_time);end
result.off_time_theor = total_off_time_theor;

if params.simulated_on
    result.burst_frequency = 1./(mean(total_on_time) + mean(total_off_time));
end
result.burst_frequency_theor = burst_frequency_theor;

if params.simulated_on && params.compute_TR;result.traveling_ratio = mean(traveling_ratio_avg);end
result.traveling_ratio_theor = traveling_ratio_theor;
if params.simulated_on
    result.burst_size_std = std(mRNA);
    temp_size = floor(size(total_off_time,2)/100);
    temp_ct_time = zeros(1,temp_size);
    for i = 1:temp_size
        temp_ct_time(i) =  1/mean(total_on_time((i-1)*100+1:i*100) + total_off_time((i-1)*100+1:i*100));
    end
    result.burst_frequency_std = std(temp_ct_time);
    result.burst_frequency_data = temp_ct_time;
end
%   result.burst_frequency_std = std(1./(total_on_time+ total_off_time);
if params.simulated_on && params.compute_TR
    result.gene_body = mean(pol_at_elong_avg);
    result.promoter = mean(pol_at_pro_avg);
    
    result.traveling_ratio_std = std(traveling_ratio_avg);
    result.promoter_std = std(pol_at_pro_avg);
    result.gene_body_std = std(pol_at_elong_avg);
end

end





