%% different parameter power-law Figure S15G-I
% ini
burst_size = [0.943628764333941 2.002311445937848 1.547553629927085 1.086367414796343];
burst_size_theor = [0.886256514774310 2.004062853184099 1.528546559304710 1.051792106623622];
burst_frequency = [0.001037326164673 0.002597802440212 0.002103798228902 0.001338097172937];
burst_frequency_std = [0.001023932805253 0.002610140613535 0.002022648360500 0.001293211376007];


x = [burst_size burst_size_theor];
y = [burst_frequency burst_frequency_std];
x = log10(x);
y = log10(y);
scatter(x,y)
axis([-0.2 0.5 -3.1 -2.5])
box on
set(gca,'TickLength',[0.02,0.025]);
t = -0.3:0.1:0.7;
s = 1.1469.*t-2.9296;
hold on
plot(t,s)
axis square

%rec
burst_size = [0.478101142299524 2.011181702668361 1.363439238238899 0.726713194482607];
burst_size_theor = [0.486483927961888 1.984686927010670 1.348905086199117 0.709355000487228];
burst_frequency = [0.001936986614364 0.002624213593928 0.002357683003437 0.002023287931758];
burst_frequency_std = [0.001934204669110 0.002636548866068 0.002368375800480 0.002052232821326];


x = [burst_size burst_size_theor];
y = [burst_frequency burst_frequency_std];
x = log10(x);
y = log10(y);
scatter(x,y)
axis([-0.5 0.6 -2.75 -2.55])
box on
set(gca,'TickLength',[0.02,0.025]);
t = -0.5:0.1:0.7;
s = 0.2219.*t-2.6450;
hold on
plot(t,s)
axis square

% rel
burst_size = [0.699778674651738 1.981271995977878 1.488771912879703 0.950045031522065];
burst_size_theor = [0.687200851406397 1.998829566237564 1.469243272507358 0.902985101534587];
burst_frequency = [0.001707690359676 0.002653308211272 0.002303010437482 0.001851414989694];
burst_frequency_std = [0.001721654157855 0.002628050208079 0.002280120268572 0.001874312038383];


x = [burst_size burst_size_theor];
y = [burst_frequency burst_frequency_std];
x = log10(x);
y = log10(y);
scatter(x,y)
 axis([-0.2 0.5 -2.8 -2.55])
box on
set(gca,'TickLength',[0.02,0.025]);
t = -0.5:0.1:0.7;
s = 0.3961.*t-2.6995;
hold on
plot(t,s)
axis square