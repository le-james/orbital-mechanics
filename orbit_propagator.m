function [state_data,time_data] = orbit_propagator(x0,delta_x0,mu,num_period,h,options)

% orbit_propagator Solves and plots the solution of the two body problem
%
%  Inputs: 
%               x: 6 states of the smaller body in km and km/s
%                  Column or row vector
%        delta_x0: initial state deviations
%              mu: standard graviational parameter of primary body
%      num_period: if custom=="default": integer number orbits to solve for
%                  if custom=="set": final time of integration
%               h: integration timestep
%         options: [custom type] 
%          custom: set a custom final time for integration: "default" or "set"  
%            type: type of dynamics to use: "twobody", "twobodystm"
% 
% Outputs:                
%          state_data: {st} or {dst stm} depends on "type" chosen
%           time_data: {h,T,final_sim_time,t_steps}
%                   T: time to complete one or more orbits
%                      for one orbit time = T/num_period
%      final_sim_time: how long the sim time is in seconds
%             t_steps: time=0 to final_sim_time in increments of h
%         t_steps_len: length of t_step
%
%            out0: output of state vectors - 6xnumber of integration steps
%            out1: output stm - 6x6xnumber of integration steps
%            out2: period of orbit
% 
% Created: July 15, 2022 by James Le - le_james@outlook.com
% Last Update: July 28, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    % sim setup

    % used to compute the orbital time with
    if all(delta_x0 == 0)
        % do nothing
    elseif options(2) == "twobody"
        x0 = x0+delta_x0;   % add deviation to initial state
    end

    if options(1) == "default"
        % orbital period
        % could reduce computation time w/o computing the full coes
        coes = cart2coes(x0,mu);
        T = coes.orbit_period;
    
        % alt code for computing the period
    %     vNorm = norm(x0(4:6));                  % velocity magnitude
    %     eng = vNorm^2/2 - mu/norm(x0(1:3));     % energy
    %     a = -mu/2/eng;                          % semi-major axis
    %     % check for hyperbolic orbit
    %     if a < 0 
    %         error("Semi-major axis is less than zero (Hyperbolic trajectory)" + ...
    %             ", Initial velocities are too high probably")
    %     end
    %     T = 2*pi*sqrt(a^3/mu);      % orbital period
    
        % number of orbits to simulate
        final_sim_time = T*num_period;
    
        % number of sim times i.e. discretize sim time with h steps
        % round up to account for integration - to fill the in the gap of the 
        % orbit trajectory plot
    %     t_steps = floor(linspace(0,final_sim_time,ceil(final_sim_time/h)));
        t_steps = linspace(0,final_sim_time,ceil(final_sim_time/h));
    
        % length of t_steps array
        t_steps_len = length(t_steps);

        % time data of default final time
        time_data = {h,T,final_sim_time,t_steps,t_steps_len};

    elseif options(1) == "set"

        % how long to integrate for
        final_sim_time = num_period;

        % number of sim times i.e. discretize sim time with h steps
        t_steps = linspace(0,final_sim_time,ceil(final_sim_time/h));

        T = "n/a";

        % length of t_steps array
        t_steps_len = length(t_steps);

        % time data of custom final time
        time_data = {h,T,final_sim_time,t_steps,t_steps_len};

    else
        error("Argument 'custom' incorrect");
    end

%     % plot orbit function - MOVED TO orbit_plotter.m
%     function plot_parameters(title_name,leg)
%         % plot orbit trajectory
% %         x = st(1,:);
% %         y = st(2,:);
% %         z = st(3,:);
% %         vx = xn(4);
% %         vy = xn(5);
% %         vz = xn(6);
% 
% %         % plot position
% %         plot3(x,y,z,'DisplayName',"hdklsjfl"); hold("on");
% 
%         title(title_name)
%         xlabel('x position [km]') 
%         ylabel('y position [km]')
%         zlabel('z position [km]') 
% 
%         [e1,e2,e3] = sphere(25);                % unit sphere
%         n = 6371;                               % earth radius
%         surf(e1*n,e2*n,e3*n,'FaceColor','w',"DisplayName","");
% 
%         legend(leg) % updates legend - doeesn't include surf of the Earth in legend
% 
%         % removes scientific notation on the plot
%         ax = gca; % axes handle
%         ax.XAxis.Exponent = 0;
%         ax.YAxis.Exponent = 0;
%         ax.ZAxis.Exponent = 0;
%         axis equal
% 
%         grid("on")
%     end

    % type of dynamics to plot
    if options(2) == "twobody"

        % plus one to account for initial states
        % n+1 more points than timesteps
        st = zeros(6,t_steps_len+1);    % state storage

        st(:,1) = x0;                   % store initial states

        for i = 1:t_steps_len
            % integrate two body dynamics
            st(:,i+1) = runge_kutta_4(@two_body_dynamics,st(:,i),mu,h);
        end

        % function output
        state_data = {st};

        % PLOTS MOVED TO orbit_plotter.m
%         % create a figure
%         figure
% 
%         % plot two body orbit
%         plot3(st(1,:),st(2,:),st(3,:)); hold("on")
% 
%         % plot figure parameters
%         plot_parameters("Two Body Orbit","Trajectory");

    elseif options(2) == "twobodystm"

        stm0 = (reshape(eye(6), 1, 36));    % initial state of stm
        st0 = cat(2,x0,stm0);               % combine dynamics and stm states

        % plus one to account for initial states
        % n+1 more points than timesteps
        st = zeros(42,t_steps_len+1);       % state storage
        stm = zeros(6,6,t_steps_len+1);     % STM matricies storage
        dx = zeros(6,t_steps_len+1);        % store propagated deviation
        dst = zeros(6,t_steps_len+1);       % store deviated states

        st(:,1) = st0;                      % store initial states

        % integrate two body dynamics and stm dynamics
        for i = 1:t_steps_len
            st(:,i+1) = runge_kutta_4(@two_body_dynamics_stm,st(:,i),mu,h);
        end

        % propagate deviation and get new deviated states
        for i = 1:t_steps_len+1                       % plus 1 to account for initial states
            stm(:,:,i) = reshape(st(7:42, i), 6, 6);  % store STM in a 3d matrix
            dx(:,i) = stm(:,:,i)*delta_x0';           % propagate deviation
            dst(:,i) = st(1:6,i)+dx(:,i);             % deviated orbit trajectory
        end

        % function output
        state_data = {dst stm};
%         out0 = dst;
%         out1 = stm;

         % PLOTS MOVED TO orbit_plotter.m
%         % output deviated states
%         out = dst;
% 
%         % output stm matrices
%         out1 = stm;
% 
%         % create a figure
%         figure
% 
%         % plot nominal position
%         plot3(st(1,:),st(2,:),st(3,:)); hold("on");
% 
%         % plot deviated position
%         plot3(dst(1,:),dst(2,:),dst(3,:)); hold("on");
% 
%         % plot figure parameters
%         leg = ["Nominal Orbit" "Deviation Orbit"];
%         plot_parameters("Two Body STM Orbit",leg);

    else
        error("Argument 'type' incorrect");
    end
    
end