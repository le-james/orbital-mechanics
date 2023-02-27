function xn_dot = cr3bp_dynamics(x,mass_ratio)

% Full nonlinear circular restricted three body problem ODE: xn_dot = cr3bp_dynamics(x,mass_ratio)
%
% NOTE THIS IS THE DIMENSIONLESS EQUATIONS!
%
%  Inputs: 
%           x: 6 states of the smaller body in km and km/s + 36 states for
%              the STM
%  mass_ratio: This is the mass ratio of the two primary bodies NOT the
%              standard gravitational parameter (the mass ratio is also
%              called mu)
% 
% Outputs:                
%      xn_dot: Differential two body and state transition matrix states 6x1 array   
% 
% To solve this ODE, need to use an integrator of your choosing
% 
% Created: February 22, 2022 by James Le - le_james@outlook.com
% Last Update: February 25, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    % store the two body states and state transition matrix states
    xn_dot = zeros(length(x),1); % reshape to 6x1

    % two body dynamics start
    
    % pull out states
    xPos = x(1);
    yPos = x(2);
    zPos = x(3);
    xVel = x(4);
    yVel = x(5);
    zVel = x(6);

    % primaries to spacecraft distance
    r1 = sqrt((xPos+mass_ratio)^2 + yPos^2 + zPos^2);
    r2 = sqrt((xPos-1+mass_ratio)^2 + yPos^2 + zPos^2);

    % potential partial derivatives
    Vx = xPos - (1-mass_ratio)*(xPos+mass_ratio)/r1^3 - mass_ratio*(xPos+mass_ratio-1)/r2^3;
    Vy = ( 1 - (1-mass_ratio)/r1^3 - mass_ratio/r2^3 )*yPos;
    Vz = -((1-mass_ratio)/r1^3 + mass_ratio/r2^3)*zPos;

    % store cr3bp states
    xn_dot(1) = xVel;
    xn_dot(2) = yVel;
    xn_dot(3) = zVel;
    xn_dot(4) = 2*yVel + Vx;
    xn_dot(5) = -2*xVel + Vy;
    xn_dot(6) = Vz;

end