function two_body_orbit_plotter(plot_data,plot_name)

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
    plot3( plot_data(1:end,1), plot_data(1:end,2), plot_data(1:end,3))

    title(plot_name);
    xlabel('x position [km]'); hold on
    ylabel('y position [km]');
    zlabel('z position [km]');

    [e1,e2,e3] = sphere(25);    % unit sphere
    n = 6371;                   % earth radius
    surf(e1*n,e2*n,-e3*n,'FaceColor','w',"DisplayName","");  % plot earth

    % wrap shpere with an image or earth
    img = imread("mario-earth.jpg");
    h = findobj("type","surface");
    set(h,"CData",img,"FaceColor","texturemap");

    ax = gca; % axes handle - get the current axis

    % removes scientific notation on the plot
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    ax.ZAxis.Exponent = 0;

    ax.Clipping = 'off';

    axis equal;

    grid("on");

end