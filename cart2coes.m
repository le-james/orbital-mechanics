function coes = cart2coes(x0,mu)
% cart2coes Converts ECI Cartesian Coordinates to Classical Keplerian 
%           Orbital Elements 
%
% Inputs: 
%           x: 6 states of the smaller body in km and km/s
%              [xPosition yPosition zPosition xVelocity yVelocity zVelocity]
%          mu: Standard graviational parameter
% 
% Outputs:                
%          coes: All classical elements plus more            
% 
% Created: July 15, 2022 by James Le - le_james@outlook.com
% Last Update: July 20, 2022
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    rVec = x0(1:3);     % pull out position vector
    vVec = x0(4:6);     % pull out velocity vector

    rNorm = norm(rVec); % magnitude of position vector
    vNorm = norm(vVec); % magnitude of velocity vector

    % energy
    eng = vNorm^2/2 - mu/rNorm;

    if eng == 0
        error("Energy is equal to zero - Parabolic orbit - Check velocities")
    elseif eng > 0
        error("Energy is greater than zero - Hyperbolic orbit - Check velocities")
    else
        % do nothing
    end

%     xHat = [1 0 0];   % not used
%     yHat = [0 1 0];   % not used
    zHat = [0 0 1];

    % specific angular momentum
    hVec = cross(rVec, vVec);
    hNorm = norm(hVec);

    % orbit inclination
    hHat = hVec/hNorm;
    i = acos(dot(zHat, hHat));

    % longitude of ascending node/ right ascension of ascending node (RAAN)
    nodeVec = cross(zHat, hVec);
    nodeNorm = norm(nodeVec);
    RAAN = acos(nodeVec(1)/nodeNorm);

    if nodeVec(2) < 0
        RAAN = 2*pi-RAAN;
    end

    % eccentricity vector
    eVec = cross(vVec, hVec)/mu - rVec/rNorm;
    eNorm = norm(eVec);

    % argument of periapsis
    nodeVec = cross(zHat, hHat);
    nodeNorm = norm(nodeVec);
    omega = acos(dot(nodeVec, eVec)/(nodeNorm*eNorm));

    if eVec(2) < 0
        omega = 2*pi-omega;
    end

    % not working
    % semi-major axis
%     a = (hNorm^2/mu)/(1-eNorm^2);
    a = -mu/2/eng;

%     % energy
%     vNorm = norm(vVec);
%     eng = vNorm^2/2 - mu/rNorm;

    % orbital period
    T = 2*pi*sqrt(a^3/mu);

    % mean motion
    n = sqrt(mu/a^3);

    % true anomaly 
    % f = acos((hNorm^2/μ/rNorm - 1)/eNorm)
    f = acos(dot(eVec, rVec)/(eNorm*rNorm));

    % eccentric anomaly
    E = 2*atan(sqrt((1-eNorm)/(1+eNorm)))*tan(f/2);

    % mean anomaly 
    M = E - eNorm*sin(E);

    % store in a struct
    coes.semi_major_axis = a;

    coes.eccentricity_vec = eVec;
    coes.eccentricity_mag = eNorm;

    coes.inclination_rad = i;
    coes.inclination_deg = rad2deg(i);

    coes.aug_periapsis_rad = omega;
    coes.aug_periapsis_deg = rad2deg(omega);

    coes.raan_rad = RAAN;
    coes.raan_deg = rad2deg(RAAN);

    coes.mean_anomaly_rad = M;
    coes.mean_anomaly_deg = rad2deg(M);

    coes.energy = eng;

    coes.angular_momentum_vec = hVec;
    coes.angular_momentum_mag = hNorm;

    coes.mean_motion = n;

    coes.orbit_period = T;

    coes.true_anomaly_rad = f;
    coes.true_anomaly_deg = rad2deg(f);

    coes.eccentric_anomaly_rad = E;
    coes.eccentric_anomaly_deg = rad2deg(E);

end