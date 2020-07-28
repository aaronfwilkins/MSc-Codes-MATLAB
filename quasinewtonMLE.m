function [BETAsave,variance,deviation,iter,BETAreturn,BETAmean,BETAstd] = quasinewtonMLE...
   (BETA0,BETAsynthetic,T_data,Y_data,Y_o,r,K,delta,burnt,limit,DesiredSD,alpha,...
   lambda,GuessSD,QNiter,options,disp_opt)
% This function performs a Maximum Likelihood Estimation (MLE) through MCMC
% and contains asymmetrical draws from a log-normal posterior
% distribution with Metropolis rejection criterion. MCMC contains a nested
% Quasi-Newton gradient root finding method to enhance the estimation
%
% Aaron Wilkins, 2018 (10076957)
%
%   - BETA0: Intial guess at the parameters
%   - T_data: Downsampled temporal set
%   - Y_data: Downsampled synthetic data / experimental data
%   - Y_o: Initial conditions (copper concentrations)
%   - r: growth rate of yeast
%   - K: copper carrying capacity of yeast culture
%   - delta: Jacobian finite difference step size (fine tuning)
%   - burnt: Burn-in time (throw-away iterations)
%   - limit: Upper/maximum iteration limit
%   - DesiredSD: standard deviation/variance of gaussian noise(fine tuning)
%   - alpha: Gradient step size (h_k in literature) (fine tuning)
%   - lambda: Tikhonov regularization coefficient (fine tuning)
%   - GuessSD: Log-normal sample multiplier (fine tuning)
%   - QNiter: Iteration limit for nested Quasi-Newton
%   - dY: system of equations
%   - options: options for selected ode solver (ode23tb)
%% =======================================================================|
% MCMC parameters and initialization
initialSD = DesiredSD;
standardev = [GuessSD;GuessSD;GuessSD;GuessSD;GuessSD;GuessSD];
BETAsave = zeros(max(size(BETA0)),(limit-burnt));
BETAmean = zeros(max(size(BETA0)),(limit-burnt));
BETAstd = zeros(max(size(BETA0)),(limit-burnt));
n = max(size(BETA0));
% Markov Chain Monte Carlo
% Metropolis-Hastings with assymetric log-normal sampling
% Maximum Likelihood Estimation (MLE) with nested Quasi-Newton Approach
% Monte Carlo Below!
%[x,fval] = runnested(BETA0,T_data,Y_data,Y_o,r,K,options,DesiredSD);
%fval
%x
iter = 1;
tic
for period = 1:10
     for mc_iter = 1:(limit/10)
        chi_square0 = sum(log((Y_data - odeSolve(T_data,Y_data, BETA0, Y_o,r,K,options)).^2)/(2*DesiredSD^2))+...  % ADDED THE NATURAL LOGARITHM TO OBTAIN THE LOG-LIKELIHOOD!!!!!!!!! 
            lambda*(norm(BETA0).^2); % Tikhonov Regularization term
        for i = 1:n
            temporary = exp(log(BETA0(i,:)) + (standardev(i,:).* randn));
            BETA1(i,:) = abs(temporary);
        end
        chi_square1 = sum(log((Y_data - odeSolve(T_data,Y_data, BETA1, Y_o,r,K,options)).^2)/(2*DesiredSD^2))+...
            lambda*(norm(BETA1).^2); % Tikhonov Regularization term
        ratio = exp(chi_square0-chi_square1);
        if rand <= ratio % Accept / Reject
            BETA0 = BETA1;
            chi_square0 = chi_square1;
            DesiredSD = DesiredSD - (DesiredSD/(limit));
            %[BETA0,fval] = runnested(BETA0,T_data,Y_data,Y_o,r,K,options,DesiredSD);
            % ================================================================|
            % Quasi-Newton Approach Below!
            [BETA0] = quasiNewton...
                (BETA0,QNiter,T_data,Y_data,Y_o,r,K,delta,alpha,options);
        end
        iter = iter + 1;
        if iter > burnt
            BETAsave(:,iter) = BETA0; %end+1
            % Store the mean per iteration
            BETAmean(1,iter) = mean(BETAsave(1,2:iter));
            BETAmean(2,iter) = mean(BETAsave(2,2:iter));
            BETAmean(3,iter) = mean(BETAsave(3,2:iter));
            BETAmean(4,iter) = mean(BETAsave(4,2:iter));
            BETAmean(5,iter) = mean(BETAsave(5,2:iter));
            BETAmean(6,iter) = mean(BETAsave(6,2:iter));
            % Store the standard deviation per iteration
            BETAstd(1,iter) = std(BETAsave(1,2:iter));
            BETAstd(2,iter) = std(BETAsave(2,2:iter));
            BETAstd(3,iter) = std(BETAsave(3,2:iter));
            BETAstd(4,iter) = std(BETAsave(4,2:iter));
            BETAstd(5,iter) = std(BETAsave(5,2:iter));
            BETAstd(6,iter) = std(BETAsave(6,2:iter)); 
        end
     end
     if disp_opt == true
         disp(['Period: ',num2str(period)])
         toc
         disp(['RMS Error: ',num2str(rms(Y_data-...
             odeSolve(T_data,Y_data, BETA0, Y_o,r,K,options)))]);
         disp(['___________________________________']);
     end
     
 end
%% =======================================================================|
% Compute statistical variance and deviation
variance = [var(BETAsave(1,:));var(BETAsave(2,:));var(BETAsave(3,:))...
    ;var(BETAsave(4,:));var(BETAsave(5,:));var(BETAsave(6,:))];
deviation = [std(BETAsave(1,:));std(BETAsave(2,:));std(BETAsave(3,:))...
    ;std(BETAsave(4,:));std(BETAsave(5,:));std(BETAsave(6,:))];
% Evaluate the mean
BETAreturn(1,:) = mean(BETAsave(1,:));
BETAreturn(2,:) = mean(BETAsave(2,:));
BETAreturn(3,:) = mean(BETAsave(3,:));
BETAreturn(4,:) = mean(BETAsave(4,:));
BETAreturn(5,:) = mean(BETAsave(5,:));
BETAreturn(6,:) = mean(BETAsave(6,:));
 end

















