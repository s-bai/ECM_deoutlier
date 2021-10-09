function [rho_min, rho_max, rho_old1_vector, rho_old2_vector, lower_asymptotes, upper_asymptotes, a0_MMA, a_MMA, c_MMA_base, d_MMA] = ...
        get_MMA_parameter(nelx, nely, n_constraint)

    %% This function prepares the MMA parameters

    rho_min = zeros(nely, nelx);
    rho_max = ones(nely, nelx);

    rho_old1_vector = zeros(nelx * nely, 1); % V0*ones(nelx*nely, 1)?
    rho_old2_vector = zeros(nelx * nely, 1); % V0*ones(nelx*nely, 1)?

    lower_asymptotes = zeros(nelx * nely, 1);
    upper_asymptotes = zeros(nelx * nely, 1);

    a0_MMA = 1;
    a_MMA = zeros(n_constraint, 1);
    c_MMA_base = ones(n_constraint, 1); % The value of c_MMA should be moderate
    d_MMA = ones(n_constraint, 1);

end
