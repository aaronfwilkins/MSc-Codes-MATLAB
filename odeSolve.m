   function [Y_test] = odeSolve(T_data,Y_data, X, Y_o,r,K,options,rshp) 
% This is an objective function that takes in a time series for
% experimental data and returns a numerical data set. This function is
% meant to be called within an inversion algorithm that supplies an 
% updated parameter set
%
% Aaron Wilkins, 2018 (10076957)
%
%    - T_data is experimental time
%    - Y_data is experimental data
%    - X is the parameters that should be varying
%    - Y_o are initial conditions
%    - r: growth rate parameter
%    - K: carrying capacity of yeast
%    - rshp: to resize Y_test
   if nargin <8
       rshp = true;
   end
        dY = @(t,Y) [r*Y(1)*(1-(Y(1)/K)); % dY(1)/dt
            X(4)^2*Y(5)-X(1)^2*Y(2); % dY(2)/dt
            X(2)^2*Y(5)-X(5)^2*Y(3)+Y(3)/(Y(3)+Y(4)+Y(5))*(r*Y(1)*...
            (1-(Y(1)/K))); % dY(3)/dt
            X(3)^2*Y(5)-X(6)^2*Y(4)+Y(4)/(Y(3)+Y(4)+Y(5))*(r*Y(1)*...
            (1-(Y(1)/K))); % dY(4)/dt
            X(1)^2*Y(2)+X(5)^2*Y(3)+X(6)^2*Y(4)-((X(3)^2+X(2)^2+X(4)^2)...
            *Y(5))+Y(5)/(Y(3)+Y(4)+Y(5))*(r*Y(1)*(1-(Y(1)/K)));]; %dY(5)/dt
         [t_test,Y_test] = ode15s(dY, T_data, Y_o,options);  
         if rshp
            Y_test = reshape(Y_test, numel(Y_test),1);
         end
   end