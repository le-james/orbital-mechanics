function [L1_x, L2_x, L3_x, L45_x, L4_y, L5_y] = compute_L_points(mass_ratio)

% Compute the five Lagrange points given the mass ratio of the system: [L1, L2, L3, L45_x, L4_y, L5_y] = compute_L_points(mass_ratio)
%
%  NOTE: This assumes z position and z velocity are zero i.e. Small
%  inclination
%
%  Inputs:                
%       mass_ratio: Mass ratio of the three body system mas_ratio = M2/(M1+M2) 
%                       M1 is the larger body, M2 is the smaller body
%  
%  Outputs: 
%             L1_x: x position of the L1 point
%             L2_x: x position of the L2 point   
%             L3_x: x position of the L3 point
%         
%             L45_x: x position of the L4 and L5 points
%             L4_y: y position of the L4 point
%             L5_y: y position of the L5 point  
% 
% 
% Created: February 22, 2022 by James Le - le_james@outlook.com
% Last Update: February 24, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    if mass_ratio > 0.5
        error("The mass ratio should be less than or equal to 0.5")
    end

    % cr3bp potential partial devrivative
    function Vx = fun(xPos)
        Vx = xPos - (1-mass_ratio)*(xPos+mass_ratio)/abs(xPos+mass_ratio)^3 - mass_ratio*(xPos+mass_ratio-1)/abs(xPos+mass_ratio-1)^3;
    end

    % intial guess to find the L points
    L1_guess = 1-mass_ratio-0.4;
    L2_guess = 1-mass_ratio+0.1;
    L3_guess = -mass_ratio-1;

    % solve for the lagrange points
    L1_x = fzero(@fun,L1_guess); % x position of the L1 point
    L2_x = fzero(@fun,L2_guess); % x position of the L2 point   
    L3_x = fzero(@fun,L3_guess); % x position of the L3 point

    L45_x = 0.5-mass_ratio; % x position of the L4 and L5 points
    L4_y = sqrt(3)/2;   % y position of the L4 point
    L5_y = -sqrt(3)/2;  % y position of the L5 point

end