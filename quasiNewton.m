function [BETA0] = quasiNewton...
    (BETA0,QNiter,T_data,Y_data,Y_o,r,K,delta,alpha,options)
% This function performs a quasi-Newton optimization method. This method is
% used when the Hessian matrix is too costly to compute, a requirement with
% the full Newton search method. The error between the data set and the
% model is evaluated and the Jacobian is used to assess local curvature and
% nudge the parameter estimation towards the closest minima
%
% Aaron Wilkins, 2018 (10076957)
%
%   - BETA0: parameter estimation
%   - QNiter: maximum iteration limit (generally around 15 to 20)
%   - T_data: Downsampled temporal data set
%   - Y_data: Downsampled synthetic / experimental data
%   - Y_o: Initial conditions (copper concentrations)
%   - r: growth rate of yeast
%   - K: copper carrying capacity of yeast culture
%   - delta: Jacobian finite difference step size (fine tuning)
%   - alpha: Gradient step size (h_k in literature) (fine tuning)
%% =======================================================================|
 for quasi = 1:QNiter
    [Y_test] = odeSolve(T_data,Y_data, BETA0, Y_o,r,K,options);
    err = Y_data - Y_test;
    jac = zeros(max(size(Y_test)),max(size(BETA0)));
    for i = 1:size(BETA0,1);
        BETA1 = BETA0;
        BETA1(i) = BETA1(i) + BETA1(i) * delta;
        jac(:,i) = (odeSolve(T_data,Y_data, BETA1, Y_o,r,K,options)-...
            Y_test)./(delta*BETA0(i));
    end
    % gradient = -2 * jac.' * err.';
    update = jac' * err;
    % delta_a = h_k * gradient;
    delta_a = (jac'*jac) \ update;
    %a = a + delta_a;
    BETA0 = BETA0 + alpha*delta_a;
 end
end