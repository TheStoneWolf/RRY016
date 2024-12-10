function [d, t_days, T_sat, P_rad, P_tot] = power_sat(a, m_focus, e, P_focus, P_gen, r_sat, sat_abs, sat_emi)

    G = 6.6743*10^(-11);
    t_elements = 100;
    E_range = 0.001;
    boltzmann = 5.67*10^-8;
    
    % Time to reach asteroid is half of the total orbital period
    t_halforbit = (2*pi*a^(3/2))/sqrt(G*m_focus)/2;
    
    t = linspace(1, t_halforbit, t_elements);
    t_days = t./(60.^2.*24);
    M = sqrt(G*m_focus/a^3)*t;
    
    E = @(E_0, e, M) E_0 - (E_0 - e.*sin(E_0)-M)./(1-e.*cos(E_0));
    E_0 = M;
    for i = 1:length(M)
        E_1 = 0;
        % This is a bit messy because there is no do-while loop in Matlab
        while E_0(i)*(1+E_range) < E_1 || E_0(i)*(1-E_range) > E_1
            E_0(i) = E_1;
            E_1 = E(E_0(i), e, M(i));
        end
    end
    d = a.*(1-e.*cos(E_0));
    J_sun_sat = P_focus ./ (4.*pi.*d.^2);
    % Assume spacecraft is a sphere with 1 m radius
    P_rad = pi .* r_sat.^2 .* sat_abs .* J_sun_sat;
    P_tot = P_rad + P_gen;
    T_sat = nthroot(P_tot ./ (4.* pi .* r_sat.^2 .* boltzmann .* sat_emi), 4); 
end

