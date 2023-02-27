function [Y0,tspan_periodic_orbit] = differential_corrector(x0_guess,mass_ratio)

% Single shooting differential corrector: [Y0,tspan_periodic_orbit] = differential_corrector(x0_guess,mass_ratio)
%  Used to find the initial conditions of a periodic orbit
%
%    Use NASA tools at https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
%    to get an intial guess for a periodic orbit
%
%    NOTE: Only Lyapunov and Halo orbits work for this differential
%    corrector
%
%  Inputs: 
%             x0_guess: Use the initial states from the NASA tools
%           mass_ratio: Mass ratio of the three body system
% 
% Outputs:                
%                   Y0: Corrected initial states of the periodic orbit and
%                   its state transition matrix - 1x42 vector
% tspan_periodic_orbit: The last element in the array is the total time of
%                       the periodic orbit in seconds
%
% Created: February 24, 2022 by James Le - le_james@outlook.com
% Last Update: February 26, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    % integrator options
    options = odeset('Events', @myEvent,"RelTol",10e-14,"InitialStep",1.0e-7,"AbsTol",10e-16);

    function [value, isterminal, direction] = myEvent(~,x)
        value = x(2);     % when y state is 0 stop the integration  
        isterminal = 1;   % Stop the integration
        direction = 0;
    end

    % tolerance for the vx and vz approx to be zero
    vel_tol = 1e-10;

    % guess periodic orbit period
    P = 5;
    
    % time steps
    tspan = 0:0.001:P;

    % initial guess with initial stm state
    Y0 = cat(2,x0_guess,reshape(eye(6),1,36));
    
    cond = 1;
    while cond
    
        % integrate ode
        [tspan_periodic_orbit,x] = ode45(@(t,x) cr3bp_dynamics_stm(x,mass_ratio),tspan,Y0,options);
    
        disp(x(end,4))
        disp(x(end,6))
    
        % check if the x and z velocities are zero
        if abs(x(end,4)) < vel_tol && abs(x(end,6)) < vel_tol
            cond = 0;   % exit loop
        end
    

        % differential correction
        
        % primaries to particle distance
        r1 = sqrt((x(end,1)+mass_ratio)^2 + x(end,2)^2 + x(end,3)^2);
        r2 = sqrt((x(end,1)-1+mass_ratio)^2 + x(end,2)^2 + x(end,3)^2);
        
        % effective potential
        Vx = x(end,1) - (1-mass_ratio)*(x(end,1)+mass_ratio)/r1^3 - mass_ratio*(x(end,1)+mass_ratio-1)/r2^3;
        Vz = -((1-mass_ratio)/r1^3 + mass_ratio/r2^3)*x(end,3);
        
        % x and z acceleration
        xdotdot = 2*x(end,5) + Vx;
        zdotdot = Vz;
        
        % stm at time T/2
        stm_T2 = reshape(x(end,7:42),6,6);
    
        % correction term
        % delta z0 and delta y0 dot correction term
        corr_terms = linsolve([stm_T2(4,3) stm_T2(4,5); stm_T2(6,3) stm_T2(6,5)] - ... 
                              [xdotdot;zdotdot]*[stm_T2(2,3) stm_T2(2,5)]/x(end,5),[-x(end,4);-x(end,6)]);
    
        % update the intial guess and repeat until x and z velocities are zero
        x0_guess = x0_guess + [0 0 corr_terms(1) 0 corr_terms(2) 0];
        
        % corrected inital guess of the periodic orbit
        Y0 = cat(2,x0_guess,reshape(eye(6),1,36));

    end
end