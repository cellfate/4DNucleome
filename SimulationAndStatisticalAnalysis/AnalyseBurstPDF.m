function TheroPDF = AnalyseBurstPDF(params)
if params.EP_flag == false
    kon1 = params.k_on1; krec = params.k_recruitment; krel = params.k_release;
    kon2 = params.k_on2; koff2 = params.k_off2;koff1 = params.k_off1;koff3 = params.k_off3;
    % BS
    Mmax = 50;
    m = 0:1:Mmax;
    theta = (krec*koff3+(krel+koff3)*koff1)/((krec+koff1)*(krel+koff3));
    BSPDF = theta.*(1-theta).^(m);

    t = 0:0.01:15000;
    % ONPDF
    x1 = 0.5.*(-(krec+krel+koff1+koff3)+sqrt((krec-krel+koff1-koff3).^2 +4.*krec.*krel));
    x2 = 0.5.*(-(krec+krel+koff1+koff3)-sqrt((krec-krel+koff1-koff3).^2 +4.*krec.*krel));
    A1 = (krec+krel+koff3+x1).*x1;
    A2 = (krec+krel+koff3+x2).*x2;
    PONPDF = -(A1.*exp(x1.*t)-A2.*exp(x2.*t))./(x1-x2);

    % OFFPDF
    y1 = 0.5.*(-(kon2+kon1+koff2)+sqrt((koff2+kon1-kon2).^2+4.*koff2.*kon2));
    y2 = 0.5.*(-(kon2+kon1+koff2)-sqrt((koff2+kon1-kon2).^2+4.*koff2.*kon2));
    B1 = (koff2+kon2+y1).*y1;
    B2 = (koff2+kon2+y2).*y2;
    POFFPDF = -(B1.*exp(y1.*t)-B2.*exp(y2.*t))./(y1-y2);

    % CTPDF
    PCTPDF = A1.*B1.*(exp(x1.*t)-exp(y1.*t))./((x1-y1).*(y1-y2).*(x1-x2))...
        +A2.*B2.*(exp(x2.*t)-exp(y2.*t))./((x2-y2).*(y1-y2).*(x1-x2))...
        -A2.*B1.*(exp(x2.*t)-exp(y1.*t))./((x2-y1).*(y1-y2).*(x1-x2))...
        -A1.*B2.*(exp(x1.*t)-exp(y2.*t))./((x1-y2).*(y1-y2).*(x1-x2));

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
    % BSPDF
    Mmax = 50;
    m = 0:1:Mmax;
    theta =  (off3.*recruitment + (release+off3).*off1)./((recruitment+off1).*(release+off3));
    BSPDFV1 = theta.*(1-theta).^(m);

    t = 0:0.01:10000;
    % ONPDF
    x1 = 0.5.*(-(recruitment+release+off1+off3)+sqrt((recruitment-release+off1-off3).^2 ...
        +4.*recruitment.*release));
    x2 = 0.5.*(-(recruitment+release+off1+off3)-sqrt((recruitment-release+off1-off3).^2 ...
        +4.*recruitment.*release));
    A1 = (recruitment+release+off3+x1).*x1;
    A2 = (recruitment+release+off3+x2).*x2;
    PONPDFV1 = -(A1.*exp(x1.*t)-A2.*exp(x2.*t))./(x1-x2);

    % OFFPDF
    y1 = 0.5.*(-(on2+on1+off2)+sqrt((off2+on1-on2).^2+4.*off2.*on2));
    y2 = 0.5.*(-(on2+on1+off2)-sqrt((off2+on1-on2).^2+4.*off2.*on2));
    B1 = (off2+on2+y1).*y1;
    B2 = (off2+on2+y2).*y2;
    POFFPDFV1 = -(B1.*exp(y1.*t)-B2.*exp(y2.*t))./(y1-y2);

    % CTPDF
    PCTPDFV1 = A1.*B1.*(exp(x1.*t)-exp(y1.*t))./((x1-y1).*(y1-y2).*(x1-x2))...
        +A2.*B2.*(exp(x2.*t)-exp(y2.*t))./((x2-y2).*(y1-y2).*(x1-x2))...
        -A2.*B1.*(exp(x2.*t)-exp(y1.*t))./((x2-y1).*(y1-y2).*(x1-x2))...
        -A1.*B2.*(exp(x1.*t)-exp(y2.*t))./((x1-y2).*(y1-y2).*(x1-x2));
    
    %% V2
    Mmax = 50;
    m = 0:1:Mmax;
    theta = (lambda_off3.*lambda_recruitment + (lambda_release+lambda_off3).*...
        lambda_off1)./((lambda_recruitment+lambda_off1).*(lambda_release+lambda_off3));
    BS_prob = zeros(size(m,2),size(d_EP,2));
    for j = 1:size(m,2)
        BS_prob(j,:) = theta.*(1-theta).^(j-1);
    end
    BSPDFV2 = (BS_prob*(P.*dt)')';

    t = 0:0.01:10000;
    % ONPDF
    x1 = 0.5.*(-(lambda_recruitment+lambda_release+lambda_off1+lambda_off3)...
        +sqrt((lambda_recruitment-lambda_release+lambda_off1-lambda_off3).^2 ...
        +4.*lambda_recruitment.*lambda_release));
    x2 =  0.5.*(-(lambda_recruitment+lambda_release+lambda_off1+lambda_off3)...
        - sqrt((lambda_recruitment-lambda_release+lambda_off1-lambda_off3).^2 ...
        +4.*lambda_recruitment.*lambda_release));
    U_ON = (lambda_recruitment+lambda_release+lambda_off3+x1).*x1;
    V_ON = (lambda_recruitment+lambda_release+lambda_off3+x2).*x2;
    
    ON_prob = zeros(size(t,2),size(d_EP,2));
    for j = 1:size(t,2)
        ON_prob(j,:) = -(U_ON.*exp(x1.*t(j))-V_ON.*exp(x2.*t(j)))./(x1-x2);
    end
    PONPDFV2 = (ON_prob*(P.*dt)')';

    % OFFPDF
    y1 = 0.5.*(-(lambda_on2+lambda_on1+lambda_off2)+...
        sqrt((lambda_off2+lambda_on1-lambda_on2).^2+4.*lambda_on2.*lambda_off2));
    y2 = 0.5.*(-(lambda_on2+lambda_on1+lambda_off2)-...
        sqrt((lambda_off2+lambda_on1-lambda_on2).^2+4.*lambda_on2.*lambda_off2));
    U_OFF = (lambda_on2+lambda_off2+y1).*y1;
    V_OFF = (lambda_on2+lambda_off2+y2).*y2;
    OFF_prob = zeros(size(t,2),size(d_EP,2));
    for j = 1:size(t,2)
        OFF_prob(j,:) = -(U_OFF.*exp(y1.*t(j))-V_OFF.*exp(y2.*t(j)))./(y1-y2);
    end
    POFFPDFV2 = (OFF_prob*(P.*dt)')';

    % CTPDF
    CT_prob = zeros(size(t,2),size(d_EP,2));
    for j = 1:size(t,2)
        CT_prob(j,:) = U_ON.*U_OFF.*(exp(x1.*t(j))-exp(y1.*t(j)))./((x1-y1).*(y1-y2).*(x1-x2))...
            +V_ON.*V_OFF.*(exp(x2.*t(j))-exp(y2.*t(j)))./((x2-y2).*(y1-y2).*(x1-x2))...
            -V_ON.*U_OFF.*(exp(x2.*t(j))-exp(y1.*t(j)))./((x2-y1).*(y1-y2).*(x1-x2))...
            -U_ON.*V_OFF.*(exp(x1.*t(j))-exp(y2.*t(j)))./((x1-y2).*(y1-y2).*(x1-x2));
    end
    PCTPDFV2 = (CT_prob*(P.*dt)')';
    
    %% mu
    minlambda = min([params.k_on1,params.k_recruitment,params.k_release]);
    maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
        params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
        99,1, 'first'))/params.friction_coef;
    mu = minlambda/maxvelocity;
    BSPDF = 1./(1+mu).*BSPDFV1+mu./(1+mu).*BSPDFV2;
    PONPDF = 1./(1+mu).*PONPDFV1+mu./(1+mu).*PONPDFV2;
    POFFPDF = 1./(1+mu).*POFFPDFV1+mu./(1+mu).*POFFPDFV2;
    PCTPDF = 1./(1+mu).*PCTPDFV1+mu./(1+mu).*PCTPDFV2;
end

if params.EP_flag == true
    TheroPDF.EP_Bin = d_EP;
    TheroPDF.EPPDF = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
end
TheroPDF.BS_Bin = m;
TheroPDF.BSPDF = BSPDF;
if params.EP_flag == true
    TheroPDF.BSPDFV1 = BSPDFV1;
    TheroPDF.BSPDFV2 = BSPDFV2;
end
TheroPDF.t = t;
TheroPDF.PONPDF = PONPDF;
TheroPDF.POFFPDF = POFFPDF;
TheroPDF.PCTPDF = PCTPDF;
if params.EP_flag == true
    TheroPDF.PONPDFV1 = PONPDFV1;
    TheroPDF.POFFPDFV1 = POFFPDFV1;
    TheroPDF.PCTPDFV1 = PCTPDFV1;
    TheroPDF.PONPDFV2 = PONPDFV2;
    TheroPDF.POFFPDFV2 = POFFPDFV2;
    TheroPDF.PCTPDFV2 = PCTPDFV2;
end

end