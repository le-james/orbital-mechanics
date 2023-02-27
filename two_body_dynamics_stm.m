function xn_dot = two_body_dynamics_stm(x,mu)

% Full nonlinear twobody dynamics state transition matrix ODE: xn_dot = two_body_dynamics_stm(~,x,mu)
%  Computes a two body trajectory and its STM
%
%  Inputs: 
%           x: 6 states of the smaller body in km and km/s + 36 states for
%              the STM
%              [xPosition yPosition zPosition xVelocity yVelocity zVelocity]
%          mu: Standard gravitational parameter of larger body
% 
% Outputs:                
%      xn_dot: Differential two body and state transition matrix states 42x1 array   
% 
% To solve this ODE, need to use an integrator of your choosing
%
% Created: July 15, 2022 by James Le - le_james@outlook.com
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

    % norm of position
    rNorm = sqrt(xPos^2 + yPos^2 + zPos^2);

    % stores next states
    xn_dot(1:3) = [xVel yVel zVel]; % velocity vector
    xn_dot(4:6) = [-mu*xPos/rNorm^3 -mu*yPos/rNorm^3 -mu*zPos/rNorm^3]; % acceleration vector

    % two body dynamics end


    
    % computing the state transition matrix, phiDot

    % jacobian of the two body ode
    G2 = zeros(6, 6);

    % top right 3x3 block of G2
    G2(1, 4) = 1;
    G2(2, 5) = 1;
    G2(3, 6) = 1;

    % low left 3x3 of G2
    G2(4, 1) = (3*mu*x(1)^2)/(rNorm^5) - mu/(rNorm^3);
    G2(5, 1) = 3*mu*x(1)*x(2)/(rNorm^5);
    G2(6, 1) = 3*mu*x(1)*x(3)/(rNorm^5);

    G2(4, 2) = 3*mu*x(1)*x(2)/(rNorm^5);
    G2(5, 2) = (3*mu*x(2)^2)/(rNorm^5) - mu/(rNorm)^(3);
    G2(6, 2) = 3*mu*x(2)*x(3)/(rNorm^5);

    G2(4, 3) = 3*mu*x(1)*x(3)/(rNorm^5);
    G2(5, 3) = 3*mu*x(2)*x(3)/(rNorm^5);
    G2(6, 3) = (3*mu*x(3)^2)/(rNorm^5) - mu/(rNorm^3);

    % form the phi states into a matrix
    phiSt = reshape(x(7:42),6,6);
    
    % propagate stm
    phiDot = G2*phiSt;

    % convert the stm from a 6x6 to a 1x36 vector
    phiDotRow = reshape(phiDot,1,36);

    % store differential stm
    xn_dot(7:42) = phiDotRow;

end