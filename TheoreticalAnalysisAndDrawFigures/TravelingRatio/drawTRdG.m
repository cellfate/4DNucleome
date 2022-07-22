%% traveling ratio for dG Figure S18
enhancer_idx = [51,55:5:95];
promoter_idx = [49,45:-5:5];
output_size = size(promoter_idx,2);

promoter =  zeros(1,output_size);
promoter_std = zeros(1,output_size);
gene_body =  zeros(1,output_size);
gene_body_std = zeros(1,output_size);
traveling_ratio = zeros(1,output_size);
traveling_ratio_std = zeros(1,output_size);
promoter_theor = zeros(1,output_size);
gene_body_theor = zeros(1,output_size);
traveling_ratio_theor = zeros(1,output_size);

addpath(genpath('.\ResultsdG'));
%% loading data
for i = 1:1:output_size
    filename = sprintf('ResultsdG\\E_%d_P_%d',enhancer_idx(i),promoter_idx(i));
    load(filename);
    
    pol_at_pro = double(results.snap_transcription_state == 4);
    pol_at_elong = results.pol_at_elong;
    pol_at_pro_all = [];
    pol_at_elong_all = [];
    for idx = 1:15
        pol_at_pro_all = [pol_at_pro_all; pol_at_pro(:,(idx-1)*666+1:idx*666)];
        pol_at_elong_all = [pol_at_elong_all; pol_at_elong(:,(idx-1)*666+1:idx*666)];
    end
    if ~sum(sum(pol_at_pro_all) == 0) && ~sum(sum(pol_at_elong_all) == 0)
        pol_at_pro_avg = sum(pol_at_pro_all)./size(pol_at_pro_all,1);
        pol_at_elong_avg = sum(pol_at_elong_all)./size(pol_at_pro_all,1);
        traveling_ratio_avg = pol_at_elong_avg./pol_at_pro_avg;
        promoter(1,i) =  mean(pol_at_pro_avg);
        promoter_std(1,i) = std(pol_at_pro_avg);
        gene_body(1,i) =   mean(pol_at_elong_avg);
        gene_body_std(1,i) = std(pol_at_elong_avg);
        traveling_ratio(1,i) = gene_body(1,i)/promoter(1,i);
        traveling_ratio_std(1,i) = std(traveling_ratio_avg);
    end
    
    
    %% lilun
    elongation_time = 100;
    burst_size = results.burst_size_theor;
    cycle_time = 1/results.burst_frequency_theor;
    gene_body_theor(1,i) = burst_size/cycle_time*elongation_time;
    traveling_ratio_theor(1,i)  = results.traveling_ratio_theor;
    promoter_theor(1,i) = gene_body_theor(1,i)/traveling_ratio_theor(1,i) ;
end

%% traveling ratio
figure3 = figure;
set(figure3,'position',[300 400 280 190],'Name','tr');
plot(enhancer_idx - promoter_idx,(traveling_ratio_theor));
hold on
scatter(enhancer_idx - promoter_idx,(traveling_ratio))
 errorbar(enhancer_idx - promoter_idx,traveling_ratio,traveling_ratio_std,'LineWidth',1);
set(gca,'TickLength',[0.02,0.025]);




