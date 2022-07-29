function [xdot,rNorm] = two_body_dynamics(x,mu)

% two_body_dynamics Full nonlinear twobody dynamics ODE
%
%  Inputs: 
%           x: 6 states of the smaller body in km and km/s
%              [xPosition yPosition zPosition xVelocity yVelocity zVelocity]
%          mu:
% 
% Outputs:                
%         xdot: Differential state
%               6x1 array       
% rNorm: Magnitude of the position vector
% 
% To solve this ODE, need to use an integrator of your choosing
%
% Created: July 15, 2022 by James Le - le_james@outlook.com
% Last Update: July 18, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    % pull out states
    xPos = x(1);
    yPos = x(2);
    zPos = x(3);
    xVel = x(4);
    yVel = x(5);
    zVel = x(6);

    % norm of position
    rNorm = sqrt(xPos^2+yPos^2+zPos^2);

    % velocity vector
    rDot = [xVel yVel zVel];

    % acceleration vector
    vDot = [-mu*xPos/rNorm^3 -mu*yPos/rNorm^3 -mu*zPos/rNorm^3];

    % stores next states
    xdot = zeros(6,1);      % reshape to 6x1
    xdot(1:3) = rDot;            
    xdot(4:6) = vDot; 

end