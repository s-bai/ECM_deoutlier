function [Alpha_objective, DAlpha_objective] = get_Alpha_objective(...
        Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, G_points, G_weights, ...
        rho_Phys, p, FE_assembly, BC, F_DOFs_and_F0, Ui, inv_T, index_concerned, phi)
    %% Calculate the CP value of the inner loop problem
    %
    %% ---------------------------------------------- Input ----------------------------------------------
    %   Lx, Ly: The nominal structural length and width
    %   nelx, nely: Elements divisions in x and y directions
    %
    %   E_bounds: The bounds of the Young's modulus samples
    %   nu_bounds: The bounds of the Poisson's ratio samples
    %
    %   theta_E0, theta_nu0: The the center of the ECM of the normalized Young's modulus and Poisson's ratio uncertainties
    %
    %   G_points, G_weights: The Gaussian quadrature points and the corresponding weighting factors
    %
    %   rho_Phys: The physical element density vector
    %
    %   p: The penalty factor
    %
    %   FE_assembly: The finite elements assembly information
    %
    %   BC: BC.free_DOF: Free DOFs
    %
    %                   ┌                                ┐
    %                   │   F1_DOF    F1_nominal_value   │
    %   F_DOFs_and_F0 = │   F2_DOF    F2_nominal_value   │
    %                   │     ...           ...          │
    %                   └                                ┘
    %
    %
    %            The arrangement of the elements of Ui
    %     ┌─────────────────┬────────────────────────┬─────────────────────────────┐
    %     │    DOFs of Ui   │ The upper limits of Ui │   The lower limits of eta   │
    %     ├─────────────────┼────────────────────────┼─────────────────────────────┤
    %     │     U1_DOF      │       U1_upper         │        eta_U1_lower         │
    %     │     U2_DOF      │       U2_upper         │        eta_U2_lower         │
    %     └─────────────────┴────────────────────────┴─────────────────────────────┘
    %
    %   inv_T: The inverse of the transformation matrix
    %    inv_T{1}: The inverse of the transformation matrix of the uncertain length and width
    %    inv_T{2}: The inverse of the transformation matrix of the uncertain Young's modulus and the Poisson's ratio
    %    inv_T{3}: The inverses of the transformation matrices of the uncertain loads
    %
    %   index_concerned: The index of the concerned structural performance
    %
    %   phi: The standard uncertain parameters vector
    %       phi(1), phi(2): For the length and width uncertainties
    %       phi(3), phi(4): For the Young's modulus and the Poisson's ratio uncertainties
    %       phi(5), phi(6): For the loads uncertainties  
    %
    %   D_KE_D_phi_flag = 0: Suppress the output of D_KE_D_phi and DF_D_phi
    %   D_KE_D_phi_flag = 1: Output D_KE_D_phi
    %
    %
    %% ---------------------------------------------- Output ----------------------------------------------
    %   Alpha_objective: The objective function value
    %
    %   DAlpha_objective: The sensitivity of the objective function

    %% Solve for displacement
    [U, ~, DUi_D_phi] = FEA(...
        Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, ...
        G_points, G_weights, rho_Phys, p, FE_assembly, BC, F_DOFs_and_F0, Ui, inv_T, phi, 1);

    %% Calculate the CP value
    Alpha_objective = U(Ui(index_concerned, 1)) - Ui(index_concerned, 2);

    %% Calculate the sensitivity of the CP value with respect to the normalized uncertain parameters
    DAlpha_objective = DUi_D_phi(:, index_concerned);

end
