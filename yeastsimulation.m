%% MCMC Statistical Model Fitting and Parameter Estimation
% This algorithm applies an MCMC algorithm to satisfy a Maximum Likelihood
% Estimation based on a frequentist approach to statistical inference and
% model fitting. Prior distributions in the parameter space are assumed to
% be uniform in this instance and the marginal likelihood is absorbed as a
% constant of proportionality. A gradient based quasi-Newton method is
% nested within the algorithm to enhance the root-finding ability
%
% Aaron Wilkins, 2018 (10076957)
% University of Calgary, 2018

clear; clc; clear all;


%=========================================================================|
% Forward Modelling synthetic data rate parameters
%

a1 =  0.3;   % forward copper reaction rate between growth media and cytosol
a2 =  0.35;  % forward reaction rate between cytosol and mitochondria
a3 =  0.3;  % forward reaction rate between cytosol and golgi network
b1 =  0.1;  % reverse reaction rate between growth media and cytosol
b2 =  0.25;  % reverse reaction rate between cytosol and mitochondria
b3 =  0.2; % reverse reaction rate between cytosol and golgi network

% NEW SET:
% a1 =  0.3;   % forward copper reaction rate between growth media and cytosol
% a2 =  0.35;  % forward reaction rate between cytosol and mitochondria
% a3 =  0.3;  % forward reaction rate between cytosol and golgi network
% b1 =  0.1;  % reverse reaction rate between growth media and cytosol
% b2 =  0.25;  % reverse reaction rate between cytosol and mitochondria
% b3 =  0.2; % reverse reaction rate between cytosol and golgi network

% OLD SET:
% a1 =  0.35;   % forward copper reaction rate between growth media and cytosol
% a2 =  0.3;  % forward reaction rate between cytosol and mitochondria
% a3 =  0.4;  % forward reaction rate between cytosol and golgi network
% b1 =  0.2;  % reverse reaction rate between growth media and cytosol
% b2 =  0.25;  % reverse reaction rate between cytosol and mitochondria
% b3 =  0.5; % reverse reaction rate between cytosol and golgi network

BETAsynthetic = [a1;a2;a3;b1;b2;b3]; % for forward modelling

%=========================================================================|
% Inverse Modelling rate parameter initial guess
%
parameterguess = 2.0;
a1 =  parameterguess;% forward copper reaction rate between growth media and cytosol
a2 =  parameterguess;% forward reaction rate between cytosol and mitochondria
a3 =  parameterguess;% forward reaction rate between cytosol and golgi network
b1 =  parameterguess;% reverse reaction rate between growth media and cytosol
b2 =  parameterguess;% reverse reaction rate between cytosol and mitochondria
b3 =  parameterguess;% reverse reaction rate between cytosol and golgi network

BETAguess = [a1;a2;a3;b1;b2;b3]; % for inversion problem

%=========================================================================|
% Define initial conditions for forward and inverse modelling
%
r = 0.2; % growth rate of yeast (also change in odesys.m)
K = 200.0; % carrying capacity of yeast (also change in odesys.m)
Y_o = [1.0; 42.0; 3.0; 5.0; 7.0]; % Initial concentration conditions
tspan = [0:0.1:60]; % time series length
sample_size = 10;
%desired_sampling_rate = 5; % Desired length of downsampled vector
%sample_size = desired_sampling_rate * 2;
DesiredSD = 0.001 ; % the desired standard deviation

%=========================================================================|
% System of ODEs
%
%=========================================================================|
% Options for selected ODE solver (ode23tb)
%
options = odeset('RelTol',1e-3,'AbsTol',1e-3);%,'MaxStep',...
     %0.5e-1,'Refine',4);
%=========================================================================|
% Forward Modelling:
%   - generate downsampled synthetic data given the conditions above
%
t_vector = [1,8,9,16,33,39,45,49,54,60];
%t_vector = [1,12,45,50,60];
%t_vector = [1,10,20,30,40,50,60];




[t2,Y2,Y2_noise,T_data,Y_data,Y_plot] = updated_datasynth(tspan,Y_o,BETAsynthetic...
    ,r,K,DesiredSD,sample_size,options,t_vector);
%[t2,Y2,Y2_noise,T_data,Y_data,Y_plot] = datasynth(tspan,Y_o,BETAsynthetic...
%    ,r,K,DesiredSD,sample_size,options);

%[x,fval] = runnested(BETAguess,T_data,Y_data,Y_o,r,K,options,DesiredSD)
figure;
plot(t2,Y2_noise(:,1),'-m','LineWidth',2); hold on; 
scatter(T_data,Y_plot(:,1),'m'); hold on; 
plot(t2,Y2_noise(:,2),'-b','LineWidth',2); hold on;
scatter(T_data,Y_plot(:,2),'b'); hold on; 
plot(t2,Y2_noise(:,3),'Color',[0.9290 0.6940 0.1250],'LineWidth',2); hold on;
scatter(T_data,Y_plot(:,3),'y'); hold on; 
plot(t2,Y2_noise(:,4),'-g','LineWidth',2); hold on;
scatter(T_data,Y_plot(:,4),'g'); hold on; 
plot(t2,Y2_noise(:,5),'-r','LineWidth',2); hold on;
scatter(T_data,Y_plot(:,5),'r'); hold on; 
ax = gca;
ax.FontSize = 24; 
title('Numerical Solutions to the System of Equations','FontSize',36);
xlabel('Time (hrs)','FontSize',36);ylabel('Copper Concentration (ng/g)','FontSize',36);
%legend({'Y(t)','N_{gro}','N_{mit}','N_{gol}','N_{cyt}'},'FontSize',30);
axis([0 60 0 210]);
%=========================================================================|
%% Inverse Modelling
%   - MCMC Quasi-Newton MLE control parameters with descriptions
%
burnt = 1; % Burn-in time (throw-away iterations)
limit = 100;% Total defined iteration limit


DesiredSD = 80.0; % Chi-square error (50.0 seemed like a good bet but wrong... in literature this value should sit at 1)
delta = 0.000001; % Jacobian finite difference step size (generally 0.0000001)
alpha = 0.1; % Gradient step size (h_k in literature)
lambda = 5; % Tikhonov regularization coefficient
GuessSD = 0.1; % Log-normal sample multiplier (generally 0.1)
QNiter = 8; % Iteration limit for nested Quasi-Newton

%=========================================================================|
% Inverse Modelling
%   - MCMC Quasi-Newton MLE function to return parameter space and stats
%
tic
[BETAsave,variance,deviation,iterations,BETAreturn,BETAmean,BETAstd] = quasinewtonMLE...
    (BETAguess,BETAsynthetic,T_data,Y_data,Y_o,r,K,delta,burnt,limit...
    ,DesiredSD,alpha,lambda,GuessSD,QNiter,options,true);
toc
time = toc;
itertime = (limit-burnt)/ time;
stats = ['(',num2str(itertime),' iterations per second)'];
disp(stats);



%% =========================================================================|
% BETA0(1,:) = mode(BETAsave(1,:));
% BETA0(2,:) = mode(BETAsave(2,:));
% BETA0(3,:) = mode(BETAsave(3,:));
% BETA0(4,:) = mode(BETAsave(4,:));
% BETA0(5,:) = mode(BETAsave(5,:));
% BETA0(6,:) = mode(BETAsave(6,:));
BETA0(1,:) = mean(BETAsave(1,0.75*limit:end));
BETA0(2,:) = mean(BETAsave(2,0.75*limit:end));
BETA0(3,:) = mean(BETAsave(3,0.75*limit:end));
BETA0(4,:) = mean(BETAsave(4,0.75*limit:end));
BETA0(5,:) = mean(BETAsave(5,0.75*limit:end));
BETA0(6,:) = mean(BETAsave(6,0.75*limit:end));
BETA0
BETAsynthetic



BETAsave(:,1)=BETAguess;
BETAsynth = BETAsave;
BETAsynth(1,:) = BETAsynthetic(1);
BETAsynth(2,:) = BETAsynthetic(2);
BETAsynth(3,:) = BETAsynthetic(3);
BETAsynth(4,:) = BETAsynthetic(4);
BETAsynth(5,:) = BETAsynthetic(5);
BETAsynth(6,:) = BETAsynthetic(6);

figure;
subplot(2,2,1);
plot(real(BETAmean(1,2:end).'),'r'); hold on;
plot(real(BETAmean(2,2:end).'),'Color',[0.9290 0.6940 0.1250]); hold on;
plot(real(BETAmean(3,2:end).'),'Color',[0.4660 0.6740 0.1880]); hold on;
plot(real(BETAmean(4,2:end).'),'r--'); hold on;
plot(real(BETAmean(5,2:end).'),'--','Color',[0.9290 0.6940 0.1250]); hold on;
plot(real(BETAmean(6,2:end).'),'--','Color',[0.4660 0.6740 0.1880]); hold on;
ax = gca;
ax.FontSize = 24;
title('Expected value per Step','FontSize',20);
xlabel('Iterations','FontSize',20);
ylabel('Expected Value (ng/g/hr)','FontSize',20);
%legend('a1','a2','a3','b1','b2','b3');
axis([0 (max(size(BETAsave))-1) min(min(BETAmean)) max(max(BETAmean))]);

subplot(2,2,2);
plot(real(BETAstd(1,2:end).'),'r'); hold on;
plot(real(BETAstd(2,2:end).'),'Color',[0.9290 0.6940 0.1250]); hold on;
plot(real(BETAstd(3,2:end).'),'Color',[0.4660 0.6740 0.1880]); hold on;
plot(real(BETAstd(4,2:end).'),'r--'); hold on;
plot(real(BETAstd(5,2:end).'),'--','Color',[0.9290 0.6940 0.1250]); hold on;
plot(real(BETAstd(6,2:end).'),'--','Color',[0.4660 0.6740 0.1880]); hold on;
ax = gca;
ax.FontSize = 24;
title('Standard Deviation per Step','FontSize',20);
xlabel('Iterations','FontSize',20);
ylabel('Standard Deviation (ng/g/hr)','FontSize',20);
%legend('a1','a2','a3','b1','b2','b3');
axis([0 (max(size(BETAsave))-1) min(min(BETAstd)) max(max(BETAstd))]);

subplot(2,2,[3,4]);
plot(real(BETAsave(1,2:end).'),'r'); hold on;
plot(real(BETAsave(2,2:end).'),'Color',[0.9290 0.6940 0.1250]); hold on;
plot(real(BETAsave(3,2:end).'),'Color',[0.4660 0.6740 0.1880]); hold on;
plot(real(BETAsave(4,2:end).'),'r--'); hold on;
plot(real(BETAsave(5,2:end).'),'--','Color',[0.9290 0.6940 0.1250]); hold on;
plot(real(BETAsave(6,2:end).'),'--','Color',[0.4660 0.6740 0.1880]); hold on;
% plot(real(BETAsynth(1,:).'),'k','LineWidth',2); hold on; %
% plot(real(BETAsynth(2,:).'),'k','LineWidth',2); hold on; %
% plot(real(BETAsynth(3,:).'),'k','LineWidth',2); hold on; %
% plot(real(BETAsynth(4,:).'),'k','LineWidth',2); hold on; %
% plot(real(BETAsynth(5,:).'),'k','LineWidth',2); hold on; %
% plot(real(BETAsynth(6,:).'),'k','LineWidth',2); hold on; %
ax = gca;
ax.FontSize = 24;
title('Parametric Convergence','FontSize',20);
xlabel('Iterations','FontSize',20);
ylabel('Parameter Values (ng/g/hr)','FontSize',20);
%legend('a1','a2','a3','b1','b2','b3');
axis([0 (max(size(BETAsave))-1) min(min(BETAsave)) max(max(BETAsave))]);





%%  Plot histograms of the rate parameters
BETAtemporary = BETAsave; % Stores BETAsave in a temporary array
BETAsave = real(BETAsave); % Takes only the real component of BETAsave
limit_plot = limit;
xlow = 0.15;
xhi = 0.8;
ylow = 0;
yhi = 2250;
wide = 0.15
% Often times, parameters corresponding to the more stiff equations do not
% find a minima until many thousand iterations have passed. In order to
% eliminate the problem of plotting irrelevant searches, a burn-in may be
% defined here and so all iterations are saved from the inversion algorithm
burnt_plot = 500;
titles = ["Forward Rate: I_{cyt}","Forward Rate: I_{mit}","Forward Rate: I_{gol}","Reverse Rate: O_{cyt}","Reverse Rate: O_{mit}","Reverse Rate: O_{gol}"];
figure;
for para = 1:min(size(BETAsave))
    BETApara = BETAsave(para,burnt_plot:limit_plot);
    stdA1 = std(BETApara);
    meanA1 = mean(BETApara);
    subplot(2,3,para)
    histogram(BETApara); hold on;
    truemeanA1 = log(mean(BETApara)/(sqrt(1+(var(BETApara)./(mean(BETApara)).^2))));
    truestdA1 = sqrt(log(1+(var(BETApara)./mean(BETApara.^2))));
    XA1 = [min(BETApara):0.0001:max(BETApara)];
    normYA1 = 1./(XA1.*truestdA1.*sqrt(2*pi)).*exp(-((log(XA1)-truemeanA1).^2)./(2*(truestdA1.^2)));
    YA1 = (max(histcounts(BETApara))/(max(normYA1)))./(XA1.*truestdA1.*sqrt(2*pi)).*exp(-((log(XA1)-truemeanA1).^2)./(2*(truestdA1.^2)));
    plot(XA1,YA1,"LineWidth",2);
    title(titles(para),'FontSize',20)
    axis([meanA1-0.5*wide meanA1+0.5*wide 0 400]);
end



















