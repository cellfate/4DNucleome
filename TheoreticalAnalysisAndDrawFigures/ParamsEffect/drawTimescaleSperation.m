% Figure S12
K = linspace(0.01,1,50);
k = linspace((0.0001),(0.01),50);
mu = zeros(length(K),length(k));
for idx = 1:1:length(K)
    for jdx = 1:1:length(k)
        input_options.EP_flag           = true;
        input_options.attraction_coef   = K(idx); % K, E-P communication coefficient
        input_options.enhancer_index    = 25; % enhancer index number
        input_options.promoter_index    = 75; % promoter index number
        
        params = ParametersBurst(input_options);
        Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
            (abs(params.promoter_index - params.enhancer_index))+params.attraction_coef)^(-1));
        dt = 0.01;
        d_EP = 0.01:dt:10;
        P = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
        
        minlambda = k(jdx);
        maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
            params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
            99,1, 'first'))/params.friction_coef;
        mu(idx,jdx) = minlambda/maxvelocity;
    end
end

figure3 = figure;
set(figure3,'position',[300 400 280 190]);
imagesc(k,K,mu);
axis xy
axis square
box on
set(gca,'TickLength',[0.02,0.025]);

K = linspace(0.1,1000,50);
k = linspace(0.0001,(0.01),50);
mu = zeros(length(K),length(k));
for idx = 1:1:length(K)
    for jdx = 1:1:length(k)
        input_options.EP_flag           = true;
        input_options.attraction_coef   = 0.1; % K, E-P communication coefficient
        input_options.enhancer_index    = 25; % enhancer index number
        input_options.promoter_index    = 75; % promoter index number
        input_options.friction_coef  = K(idx); % friction coefficient [Newton][sec]/[mu m]
        input_options.diffusion_const = 4e-3/input_options.friction_coef; % diffusion constant [mu m]^2 /[sec]
        
        params = ParametersBurst(input_options);
        Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
            (abs(params.promoter_index - params.enhancer_index))+params.attraction_coef)^(-1));
        dt = 0.01;
        d_EP = 0.01:dt:10;
        P = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
        
        minlambda = k(jdx);
        maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
            params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
            99,1, 'first'))/params.friction_coef;
        mu(idx,jdx) = minlambda/maxvelocity;
    end
end

figure3 = figure;
set(figure3,'position',[300 400 280 190]);
imagesc(k,K,mu);
axis xy
axis square
box on
set(gca,'TickLength',[0.02,0.025]);

