%% traveling ratio for Kep Figure S18
K_coefficient1 =  [ -1.8000   -1.6000   -1.4000   -1.2000   -1.0000  -0.8000   -0.6000   -0.4000   -0.2000   -0.0000];
output_size = size(K_coefficient1,2);

traveling_ratio = zeros(1,output_size);
traveling_ratio_std = zeros(1,output_size);
promoter_theor = zeros(1,output_size);
gene_body_theor = zeros(1,output_size);
traveling_ratio_theor = zeros(1,output_size);

addpath(genpath('.\ResultsKEP'));
%% loading data
for i = 1:1:output_size
   
    filename = sprintf('ResultsKEP\\K_%d.mat',10^(K_coefficient1(i)));
    load(filename);
    pol_at_pro = double(results.snap_transcription_state == 4);
    pol_at_elong = results.pol_at_elong;
    pol_at_pro_all = [];
    pol_at_elong_all = [];
    for idx = 1:20
        pol_at_pro_all = [pol_at_pro_all; pol_at_pro(:,(idx-1)*500+1:idx*500)];
        pol_at_elong_all = [pol_at_elong_all; pol_at_elong(:,(idx-1)*500+1:idx*500)];
    end
    if ~sum(sum(pol_at_pro_all) == 0) && ~sum(sum(pol_at_elong_all) == 0)
        pol_at_pro_avg = sum(pol_at_pro_all)./size(pol_at_pro_all,1);
        pol_at_elong_avg = sum(pol_at_elong_all)./size(pol_at_pro_all,1);
        traveling_ratio_avg = pol_at_elong_avg./pol_at_pro_avg;
        traveling_ratio(1,i) = mean(traveling_ratio_avg);
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
plot(K_coefficient1,(traveling_ratio_theor));
hold on
scatter(K_coefficient1,(traveling_ratio))
 errorbar(K_coefficient1,traveling_ratio,traveling_ratio_std*0.5,'LineWidth',1);
set(gca,'TickLength',[0.02,0.025]);




