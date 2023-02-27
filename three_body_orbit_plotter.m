function three_body_orbit_plotter(plot_data,mass_ratio,plot_name,options)

% orbit_plotter Plots the states of one or more orbits
%
%  Inputs: 
%           plot_data: states of two body orbit
%           plot_name: name of the plot
%
% Outputs:                
%           Void function that outputs a plot of an orbit
% 
% Created: July 21, 2022 by James Le - le_james@outlook.com
% Last Update: February 23, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    % plot the states
    plot3(plot_data(1:end,1), plot_data(1:end,2), plot_data(1:end,3), DisplayName='Trajectory')

    title(plot_name);
    xlabel('x position'); hold on
    ylabel('y position');
    zlabel('z position');

    if any(options) == 1

        % solve for Lagrange points
        [L1_x, L2_x, L3_x, L45_x, L4_y, L5_y] = compute_L_points(mass_ratio);
        
        % plot of lagrange points and primaries
        if options(1) == 1
            plot(-mass_ratio,0,Marker="*",MarkerFaceColor="auto",LineWidth=2,MarkerSize=10,DisplayName='Primary')
        end

        if options(2) == 1
            plot(1-mass_ratio,0,Marker="+",MarkerFaceColor="auto",LineWidth=2,MarkerSize=10,DisplayName='Secondary')
        end

        if options(3) == 1
            plot(L1_x,0,Marker="o",MarkerFaceColor="auto",LineWidth=2,DisplayName='L1')
        end

        if options(4) == 1
            plot(L2_x,0,Marker="o",MarkerFaceColor="auto",LineWidth=2,DisplayName='L2')
        end

        if options(5) == 1
            plot(L3_x,0,Marker="o",MarkerFaceColor="auto",LineWidth=2,DisplayName='L3')
        end

        if options(6) == 1
            plot(L45_x,L4_y,Marker="o",MarkerFaceColor="auto",LineWidth=2,DisplayName='L4')
        end

        if options(7) == 1
            plot(L45_x,L5_y,Marker="o",MarkerFaceColor="auto",LineWidth=2,DisplayName='L5')
        end

    end

    % display legend
    legend();

    ax = gca; % axes handle - get the current axis

    % removes scientific notation on the plot
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    ax.ZAxis.Exponent = 0;

    ax.Clipping = 'off';

%     h = get(gca,'DataAspectRatio');
%     if h(3) == 1
%           set(gca,'DataAspectRatio',[2 2 2/max(2*h(1:2))])
%     else
%           set(gca,'DataAspectRatio',[2 2 2*h(3)])
%     end

    grid("on");

end