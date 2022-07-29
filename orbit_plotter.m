function orbit_plotter(plot_data)

% orbit_plotter Plots the states of one or more orbits
%
%  Inputs: 
%           plot_data: cell array containing states
%
% Outputs:                
%           Void function that outputs a plot of an orbit
%
% Example:
%           plot_data = {orbit_1_states "orbit_1_name" orbit_2_states "orbit_2_name" ...}
%           orbit_plotter(plot_data)
% 
% Created: July 21, 2022 by James Le - le_james@outlook.com
% Last Update: July 28, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    leg_name = strings(1,length(plot_data)/2);     % name of plots

    k = 0;  % every odd element of plot_data

    for i = 1:length(plot_data)/2

        orbit_pos = plot_data{i+k};
        
        % pull out position
        x = orbit_pos(1,:);
        y = orbit_pos(2,:);
        z = orbit_pos(3,:);

        % pull out velocity
    %     vx = xn(4);
    %     vy = xn(5);
    %     vz = xn(6);

%         leg_name(i) = sprintf("Orbit %d",i);
        leg_name(i) = plot_data{i+k+1};

        % plot orbit position trajectory
        plot3(x,y,z); hold("on");

        k=k+1; % every odd element of plot_data

    end

    title("Two Body Orbit(s)");
    xlabel('x position [km]');
    ylabel('y position [km]');
    zlabel('z position [km]');

    [e1,e2,e3] = sphere(25);    % unit sphere
    n = 6371;                   % earth radius
    surf(e1*n,e2*n,-e3*n,'FaceColor','w',"DisplayName","");  % plot earth

    % wrap shpere with an image or earth
    img = imread("mario-earth.jpg");
    h = findobj("type","surface");
    set(h,"CData",img,"FaceColor","texturemap");

    legend(leg_name);    % legend of plot

    % removes scientific notation on the plot
    ax = gca; % axes handle - get the current axis
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    ax.ZAxis.Exponent = 0;
    ax.Clipping = 'off';
    axis equal;

    grid("on");

end