function [E_actual, nu_actual] = get_E_nu(E_bounds, nu_bounds, theta_E0, theta_nu0, phi, inv_T)
    %% This local function calculates the actual Young's modulus and Poisson's ratio
    %
    %% ---------------------------------- Input ----------------------------------
    %   E_bounds, nu_bounds: The bounds of the Young's modulus and Poisson's ratio samples
    %
    %   theta_E0, theta_nu0: The center of the ECM of the normalized Young's modulus and Poisson's ratio uncertainties
    %
    %   phi: The standardized uncertain parameters vector
    %       phi(3), phi(4): The phi modelling the Young's modulus and the Poisson's ratio uncertainties
    %
    %   inv_T:  
    %    inv_T{2}: The inverse of the transformation matrix of the uncertain Young's modulus and the Poisson's ratio
    %
    %% ---------------------------------- Output ----------------------------------
    %   E_actual: The actual Young's modulus
    %   nu_actual: The actual Poisson's ratio

    %%
    theta_E_theta_nu = inv_T{2} * phi(3:4) + [theta_E0; theta_nu0];

    theta_E = theta_E_theta_nu(1);
    theta_nu = theta_E_theta_nu(2);

    E_actual = theta_E * (E_bounds(2) - E_bounds(1)) + E_bounds(1);
    nu_actual = theta_nu * (nu_bounds(2) - nu_bounds(1)) + nu_bounds(1);
end
