function xn_dot = cr3bp_dynamics_stm(x,mass_ratio)
    
% Full nonlinear circular restricted three body problem state transition matrix ODE: xn_dot = cr3bp_dynamics(x,mass_ratio)
%  Computes a three body trajectory and its STM in the rotating coordinate system
%
%  NOTE: THIS IS THE DIMENSIONLESS EQUATIONS!
%
%  Inputs: 
%           x: 6 states of the smaller body in km and km/s + 36 states for
%              the STM
%              [xPosition yPosition zPosition xVelocity yVelocity zVelocity]
%  mass_ratio: This is the mass ratio of the two primary bodies NOT the
%              standard gravitational parameter (the mass ratio is also
%              called mu)
% 
% Outputs:                
%      xn_dot: Differential two body and state transition matrix states 42x1 array   
% 
% To solve this ODE, need to use an integrator of your choosing
% 
% Created: February 22, 2022 by James Le - le_james@outlook.com
% Last Update: February 23, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    % store the two body states and state transition matrix states
    xn_dot = zeros(length(x),1); % reshape to 42x1

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

    % double partial derivatives to the cr3bp potential
    Vxx = (mass_ratio - 1)/((mass_ratio + xPos)^2 + yPos^2 + zPos^2)^(3/2) - mass_ratio/((mass_ratio + xPos - 1)^2 + yPos^2 + zPos^2)^(3/2) + (3*mass_ratio*(2*mass_ratio + 2*xPos - 2)^2)/(4*((mass_ratio + xPos - 1)^2 + yPos^2 + zPos^2)^(5/2)) - (3*(2*mass_ratio + 2*xPos)^2*(mass_ratio - 1))/(4*((mass_ratio + xPos)^2 + yPos^2 + zPos^2)^(5/2)) + 1;
    Vyy = (mass_ratio - 1)/((mass_ratio + xPos)^2 + yPos^2 + zPos^2)^(3/2) - mass_ratio/((mass_ratio + xPos - 1)^2 + yPos^2 + zPos^2)^(3/2) - (3*yPos^2*(mass_ratio - 1))/((mass_ratio + xPos)^2 + yPos^2 + zPos^2)^(5/2) + (3*mass_ratio*yPos^2)/((mass_ratio + xPos - 1)^2 + yPos^2 + zPos^2)^(5/2) + 1;
    Vzz = (mass_ratio - 1)/((mass_ratio + xPos)^2 + yPos^2 + zPos^2)^(3/2) - mass_ratio/((mass_ratio + xPos - 1)^2 + yPos^2 + zPos^2)^(3/2) - (3*zPos^2*(mass_ratio - 1))/((mass_ratio + xPos)^2 + yPos^2 + zPos^2)^(5/2) + (3*mass_ratio*zPos^2)/((mass_ratio + xPos - 1)^2 + yPos^2 + zPos^2)^(5/2);
    Vxy = (3*mass_ratio*yPos*(2*mass_ratio + 2*xPos - 2))/(2*((mass_ratio + xPos - 1)^2 + yPos^2 + zPos^2)^(5/2)) - (3*yPos*(2*mass_ratio + 2*xPos)*(mass_ratio - 1))/(2*((mass_ratio + xPos)^2 + yPos^2 + zPos^2)^(5/2));
    Vxz = (3*mass_ratio*zPos*(2*mass_ratio + 2*xPos - 2))/(2*((mass_ratio + xPos - 1)^2 + yPos^2 + zPos^2)^(5/2)) - (3*zPos*(2*mass_ratio + 2*xPos)*(mass_ratio - 1))/(2*((mass_ratio + xPos)^2 + yPos^2 + zPos^2)^(5/2));
    Vyz = (3*mass_ratio*yPos*zPos)/((mass_ratio + xPos - 1)^2 + yPos^2 + zPos^2)^(5/2) - (3*yPos*zPos*(mass_ratio - 1))/((mass_ratio + xPos)^2 + yPos^2 + zPos^2)^(5/2);

    % stm dynamics jacobian - G2 aka A
    G2 = [ 0   0   0   1 0 0;
           0   0   0   0 1 0;
           0   0   0   0 0 1;
          Vxx Vxy Vxz  0 2 0;
          Vxy Vyy Vyz -2 0 0;
          Vxz Vyz Vzz  0 0 0];
    
    % form the phi states into a matrix
    phiSt = reshape(x(7:42),6,6);
    
    phiDot = G2*phiSt;  % propagate stm
    phiDotRow = reshape(phiDot,1,36);
    
    % store cr3bp stm states
    xn_dot(7:42) = phiDotRow; % differential stm
end