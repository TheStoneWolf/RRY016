m_earth = 5.9722*10^24;
m_sun = 2*10^30;
G = 6.6743*10^(-11);
r_earth = 6378;
r_au = 149.597871*10^9;
r_asteroid = 2.5 * r_au;
r_sun = 696340*1000;
v_orbit_e = 29.8 * 1000;
T_sun = 5800;
P_sun = 3.86*10^26;
P_gen = 150;
sat_abs = 0.3; sat_emi = 0.6;
albedo = 0.34; F = 0.3;
r_sat = 1;
boltzmann = 5.67*10^-8;

my_fsize = 15;

v_closed = @(M, r, a) sqrt(G.*M.*(2./(1000*r)-1./(1000*a)));

% Hohmann transfer earth-asteroid directly
close all;
r_park = r_earth + 2000;
t_elements = 1000;
E_range = 0.001;

a_earth_asteroid = (r_asteroid + r_au + r_park) / 2;
v_escape = sqrt(2*G*m_earth/r_park);
v_natural = v_closed(m_earth, r_park, r_park);
v_earth_asteroid = v_closed(m_sun, r_au + r_park, a_earth_asteroid);
delta_1 = v_earth_asteroid + v_escape - v_natural

e = r_asteroid/a_earth_asteroid-1;

%% Use example P_gen and r_sat
[d_sun, t_days, T_sat, P_rad, P_tot] = power_sat(a_earth_asteroid, ...
                m_sun, e, P_sun, P_gen, r_sat, sat_abs, sat_emi);

    figure(1)
    plot(t_days, d_sun./r_au), title("Distance from sun"),
    xlabel("t (days)"), ylabel("d (AU)")

    fontsize(my_fsize, "points")
    
    figure(2)
    plot(t_days, P_rad, 'color', 'yellow'), title("Heat gen per day"),
    xlabel("t (days)"), ylabel("P (W)"), hold on
    plot(t_days, P_gen.*ones(length(t_days),1), 'color', 'blue')
    plot(t_days, P_tot, 'color', 'red')
    legend('Sun', 'Satellite', 'Total'), fontsize(16, 'points')
    fontsize(my_fsize, "points")
    hold off
    
    figure(3)
    plot(t_days, T_sat - 273.15), title("Temperature per day"),
    fontsize(my_fsize, "points")
    xlabel("t (days)"), ylabel("T (C)");

%% Plot for different P_gen and r_sat

ele = 100;
P_gens = linspace(60, 2000, ele);
r_sats = linspace(0.2, 5, ele);
%d_sun = zeros(ele, ele); t_days = zeros(ele);
T_min = zeros(ele, ele); T_sat_whole = zeros(ele);
T_max = zeros(ele, ele);
P_rad = zeros(ele); P_tot = zeros(ele);
share_P = zeros(ele, ele);

for i = 1:length(P_gens)
    for j = 1:length(r_sats)
        
        [~, ~, T_sat_whole, P_rad, P_tot] = power_sat(a_earth_asteroid, ...
            m_sun, e, P_sun, P_gens(i), r_sats(j), sat_abs, sat_emi);
        % Don't know why flipping like this is needed
        T_min(j,i) = T_sat_whole(end) -273.15;

        J_albedo = albedo*F*P_rad(1);
        J_thermal = 237*(r_earth/r_park)^2;
        A = pi*r_sats(j)^2;
        P_extra = A*(sat_abs*J_albedo+sat_emi*J_thermal)+P_tot(1);
        P_max = P_extra + P_tot(1);
        share_P(j,i) = P_extra./P_max;
        T_max(j,i) = nthroot(P_max/(sat_emi*A*boltzmann),4)-273.15;
    end
end

figure(4)
[X, Y] = meshgrid(P_gens, r_sats);
[C, h] = contour(X, Y, T_min, 'ShowText','on'); title("min T during journey"), hold on
clabel(C, h, 'FontSize', my_fsize)
fontsize(my_fsize, "points")
xlabel("Power Generation (W)")
ylabel("Spacecraft Radius (m)")
plot(P_gen, r_sat, 'r.', 'MarkerSize', 2*my_fsize)
hold off

figure(5)
[C, h] = contour(X, Y, T_max, 'ShowText', 'on'); title("T at orbit"), hold on
clabel(C, h, 'FontSize', my_fsize)
fontsize(my_fsize, "points")
xlabel("Power Generation (W)")
ylabel("Spacecraft Radius (m)")
plot(P_gen, r_sat, 'r.', 'MarkerSize', 2*my_fsize)
hold off

figure(6)
[C, h] = contour(X, Y, share_P, 'ShowText', 'on'); title("Share of P from earth"), hold on
clabel(C, h, 'FontSize', my_fsize)
fontsize(my_fsize, "points")
xlabel("Power Generation (W)")
ylabel("Spacecraft Radius (m)")
plot(P_gen, r_sat, 'r.', 'MarkerSize', 2*my_fsize)
hold off