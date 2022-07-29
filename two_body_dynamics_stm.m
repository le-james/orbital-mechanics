function xn = two_body_dynamics_stm(x,mu)

% two_body_dynamics_stm Full nonlinear twobody dynamics ODE
%
%  Inputs: 
%           x: 6 states of the smaller body in km and km/s
%              [xPosition yPosition zPosition xVelocity yVelocity zVelocity]
%          mu:
% 
% Outputs:                
%           xn: Differential state
% 
% To solve this ODE, need to use an integrator of your choosing
%
% Created: July 15, 2022 by James Le - le_james@outlook.com
% Last Update: July 18, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    % output of two body dynamics needed for stm computation
    [two_body_states,rNorm] = two_body_dynamics(x(1:6),mu);

    % stm dynamics jacobian
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
    phiSt = reshape(x(7:42), 6, 6);
    
    % get stm ode and reshape to row
%     phiDotRow = G2;                   % same results as below - 
                                        % only slightly diff numbers in the stm
    phiDot = G2*phiSt;                  % propagate stm
    phiDotRow = reshape(phiDot, 1, 36);
    
    % store states
    xn = zeros(length(two_body_states)+length(G2)^2,1); % reshape to 42x1

    xn(1:6) = two_body_states;  % two body position and velocity

    xn(7:42) = phiDotRow;       % differential stm

end