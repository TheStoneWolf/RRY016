m_earth = 5.9722*10^24;
m_sun = 2*10^30;
G = 6.6743*10^(-11);
r_earth = 6378;
r_au = 149.597871*10^6;
r_asteroid = 2.5 * r_au;
v_orbit_e = 29.8 * 1000;

v_closed = @(M, r, a) sqrt(G.*M.*(2./(1000*r)-1./(1000*a)));

%% Hohmann transfer earth-sun-asteroid
r_park = r_earth + 2000; % 2000 km above surface
r_dist_sun = 10*10^6; % 10 million km

% Calculations

v_parked = v_closed(m_earth, r_park, r_park) + v_orbit_e;
a_earth_sun = (r_dist_sun + r_au + r_park)/2;
v_earth_sun_e = v_closed(m_sun, r_park, a_earth_sun);
delta_1 = v_earth_sun_e - v_parked

v_earth_sun_s = v_closed(m_sun, r_dist_sun, a_earth_sun);
a_asteroid_sun = r_asteroid + r_dist_sun;
v_asteroid_sun = v_closed(m_sun, r_dist_sun, a_asteroid_sun);
delta_2 = v_asteroid_sun - v_earth_sun_s

%% Hohmann transfer earth-asteroid directly
r_park = r_earth + 2000;

a_earth_asteroid = (r_asteroid + r_au + r_park) / 2;
v_parked = v_closed(m_earth, r_park, r_park) + v_orbit_e
v_earth_asteroid = v_closed(m_sun, r_au + r_park, a_earth_asteroid)
delta_1 = v_earth_asteroid - v_parked

%% Hohmann transfer earth-asteroid directly (inverse)
r_park = r_earth + 2000;

a_earth_asteroid = (r_asteroid + r_au - r_park) / 2;
v_parked = -v_closed(m_earth, r_park, r_park) + v_orbit_e
v_earth_asteroid = v_closed(m_sun, r_au - r_park, a_earth_asteroid)
delta_1 = v_earth_asteroid - v_parked


