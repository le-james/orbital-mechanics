function xn = runge_kutta_4(callback_function,x,mu,h)

% runge_kutta_4 Integrates an Ordinary Differential Equation (ODE) using 
%              Runge-Kutta 4 method
%  
% Inputs: 
%           callback_function: ODE to integrate
%                           x: State to integrate from e.g. Initial
%                              condition
%                h (timestep): Unit of seconds usually
%                              It's the "resolution" of the integration
%                              The smaller the timestep the more accurate
% 
% Outputs:                
%                          xn: Integrated state after one timestep
%                              The xn array shape is determined on the
%                              output array shape of the callback function
% 
% Example:
%           x = [ ... ];    % array of inital states/conditions
%           timestep = #;   % timestep of your choosing [seconds]
%           xn = runge_kutta_4(@my_function,x,timestep);
%
% Created: July 15, 2022 by James Le - le_james@outlook.com
% Last Update: July 18, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    % runge kutta 4 algorithm
    f1 = callback_function(x,mu);
    f2 = callback_function(x + 0.5*h*f1,mu);
    f3 = callback_function(x + 0.5*h*f2,mu);
    f4 = callback_function(x + h*f3,mu);
    xn = x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4);
    
end