%% two body orbit simulatation example given an initial position and velocity
clear; clc;

addpath("../")

% earth standard gravitational parameter
mu = 398600.4418;

% number of orbit periods to sim
num_periods = 1;

% integration timestep - h = 1min = 60sec
h = 60;

% initial conditions
v_circular = sqrt(mu/10000);       % circular orbit velocity
x0 = [10000 0 0 0 v_circular 0];   % initial state km and km/s
dx0_1 = [0 0 0 0 0 0];             % no deviation to orbit
dx0_2 = [0 0 0 0 0 1];             % deviation in z velocity by 1km/s

% simulate two orbits
options = ["default" "twobody"];
[orbit1_data,time1_data] = orbit_propagator(x0,dx0_1,mu,num_periods,h,options);
[orbit2_data,time2_data] = orbit_propagator(x0,dx0_2,mu,num_periods,h,options);

% position and velocity of orbit
st1 = orbit1_data{1};
st2 = orbit2_data{1};

% plot both orbits
plot_data = {st1 "Orbit 1" st2 "Orbit 2"};
orbit_plotter(plot_data);

% initial position marker
plot3(x0(1),x0(2),x0(3),".",'MarkerSize',35,DisplayName="Initial Position")
