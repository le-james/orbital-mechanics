function xn_dot = two_body_dynamics(x,mu)

% Full nonlinear twobody dynamics ODE: xn_dot = two_body_dynamics(~,x,mu)
%
%  Inputs: 
%           x: 6 states of the smaller body in km and km/s
%              [xPosition yPosition zPosition xVelocity yVelocity zVelocity]
%          mu: Standard gravitational parameter of larger body
% 
% Outputs:                
%      xn_dot: Differential states 6x1 array   
%
% rNorm: Magnitude of the position vector
% 
% To solve this ODE, need to use an integrator of your choosing
%
% Created: July 15, 2022 by James Le - le_james@outlook.com
% Last Update: February 22, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    % pull out states
    xPos = x(1);
    yPos = x(2);
    zPos = x(3);
    xVel = x(4);
    yVel = x(5);
    zVel = x(6);

    % norm of position
    rNorm = sqrt(xPos^2 + yPos^2 + zPos^2);

    % velocity vector
    rDot = [xVel yVel zVel];

    % acceleration vector
    vDot = [-mu*xPos/rNorm^3 -mu*yPos/rNorm^3 -mu*zPos/rNorm^3];

    % stores next states
    xn_dot = zeros(6,1);      % reshape to 6x1
    xn_dot(1:3) = rDot;            
    xn_dot(4:6) = vDot; 

end