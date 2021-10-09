function [U, K, D_Ui_D_phi, D_Ui_D_rho, c, D_c_D_phi, D_c_D_rho] = FEA(...
        Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, ...
        G_points, G_weights, rho_Phys, p, FE_assembly, BC, F_DOFs_and_F0, Ui, inv_T, phi, D_KE_D_phi_flag)

    %% Finite element analysis subroutine

    %     Input arguments                                                           Description
    %  --------------------- --------------------------------------------------------------------------------------------------------------------------
    %   Lx                    The nominal length
    %   Ly                    The nominal width
    %   nelx, nely            Element numbers in x and y directions
    %   E_bounds              The bounds of the Young's modulus samples
    %   nu_bounds             The bounds of the Poisson's ratio samples
    %   theta_E0, theta_nu0   The the center of the ECM of the normalized Young's modulus and Poisson's ratio uncertainties
    %   G_points              Two-dimensional Gaussian points
    %   G_weights             Two-dimensional Gaussian weights (tensor product)
    %   rho_Phys              The physical element density vector
    %   p                     The penalty factor
    %   FE_assembly           The finite elements assembly information
    %   BC.free_DOF           Free DOFs
    %   F_DOFs_and_F0         The DOFs and the nominal values of loads (see the README for details)
    %   Ui                    The concerned displacements' DOFs, upper limits and the lower limits on reliability indices (See the README for details)
    %   inv_T{1}              The inverse of the transformation matrix of the uncertain length and width
    %   inv_T{2}              The inverse of the transformation matrix of the uncertain Young's modulus and the Poisson's ratio
    %   inv_T{3}              The inverse of the transformation matrix of the uncertain loads
    %   phi                   The standard uncertain parameters vector
    %   phi(1), phi(2)        For the length and width uncertainties
    %   phi(3), phi(4)        For the Young's modulus and the Poisson's ratio uncertainties
    %   D_KE_D_phi_flag = 0   Suppress the output of D_KE_D_phi
    %   D_KE_D_phi_flag = 1   Output D_KE_D_phi

    %   Output arguments                                           Description
    %  ------------------ ----------------------------------------------------------------------------------------------
    %   U                  The displacements
    %   D_Ui_D_phi         Sensitivity of the concerned displacements with respect to the standard uncertain parameters
    %   D_Ui_D_rho         Sensitivity of the concerned displacements with respect to the element densities
    %   K                  The stiffness matrix
    %   c                  The mean compliance
    %   D_c_D_phi          Sensitivity of the mean compliance with respect to the standard uncertain parameters
    %   D_c_D_rho          Sensitivity of the mean compliance with respect to the element densities

    %% Setup the finite element parameters
    % The lower limit of the Young's modulus
    Emin = 1e-9;

    %% Calculate the actual Young's modulus
    [E_actual, ~] = get_E_nu(E_bounds, nu_bounds, theta_E0, theta_nu0, phi, inv_T);

    %% Get element stiffness matrix KE (NOT multiplied by the Young's modulus E)
    if D_KE_D_phi_flag ~= 0
        [KE, D_KE_D_phi_xy, D_KE_D_E, D_KE_D_nu] = get_KE(...
            Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, phi, inv_T, G_points, G_weights);
    else
        [KE] = get_KE(...
            Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, phi, inv_T, G_points, G_weights);

        D_c_D_phi = 0;
    end

    %% Get the sparse load vector
    % The number of uncertain loads components
    [n_F, ~] = size(F_DOFs_and_F0);

    % The actual loads components
    F_actual = inv_T{3} * phi(5:end) + F_DOFs_and_F0(:, 2);

    % The actual sparse load vector
    F = sparse(F_DOFs_and_F0(:, 1), ones(size(F_DOFs_and_F0(:, 1))), F_actual, 2 * (nelx + 1) * (nely + 1), 1);

    %% Solve the FEA equation for the displacements
    iK = FE_assembly.iK;
    jK = FE_assembly.jK;

    sK = reshape(KE(:) * (Emin + rho_Phys(:).'.^p * (E_actual - Emin)), 64 * nelx * nely, 1);

    K = sparse(iK, jK, sK);
    K = (K + K.') / 2;

    if D_KE_D_phi_flag ~= 0
        D_sK_D_phi_E = reshape(D_KE_D_E(:) * rho_Phys(:).'.^p * (E_bounds(2) - E_bounds(1)) * inv_T{2}(1, 1), 64 * nelx * nely, 1) + ...
            reshape(D_KE_D_nu(:) * (Emin + rho_Phys(:).'.^p * (E_actual - Emin)) * (nu_bounds(2) - nu_bounds(1)) * inv_T{2}(2, 1), 64 * nelx * nely, 1);

        D_sK_D_phi_nu = reshape(D_KE_D_E(:) * rho_Phys(:).'.^p * (E_bounds(2) - E_bounds(1)) * inv_T{2}(1, 2), 64 * nelx * nely, 1) + ...
            reshape(D_KE_D_nu(:) * (Emin + rho_Phys(:).'.^p * (E_actual - Emin)) * (nu_bounds(2) - nu_bounds(1)) * inv_T{2}(2, 2), 64 * nelx * nely, 1);

        D_K_D_phi_E = sparse(iK, jK, D_sK_D_phi_E);
        D_K_D_phi_E = (D_K_D_phi_E + D_K_D_phi_E.') / 2;

        D_K_D_phi_nu = sparse(iK, jK, D_sK_D_phi_nu);
        D_K_D_phi_nu = (D_K_D_phi_nu + D_K_D_phi_nu.') / 2;

    end

    % Free DOFs
    free_DOFs = BC.free_DOFs;

    % Solve the FEA equation using the back slash operator "\"
    U = zeros(2 * (nelx + 1) * (nely + 1), 1);
    U(free_DOFs) = K(free_DOFs, free_DOFs) \ F(free_DOFs);

    %% Calculate the sensitivity of the concerned displacements with respect to the standard uncertain parameters
    % Initialize the adjoint vectors
    adjoint_lambda = zeros(2 * (nelx + 1) * (nely + 1), 2);

    % Initialize the unit vectors
    unit_vector = zeros(2 * (nelx + 1) * (nely + 1), 2);

    for ii = 1:2
        unit_vector(Ui(ii, 1), ii) = 1;

        adjoint_lambda(free_DOFs, ii) = K(free_DOFs, free_DOFs) \ unit_vector(free_DOFs, ii);
    end

    %% Initialization
    D_Ui_D_phi = zeros(4 + n_F, 2);

    if D_KE_D_phi_flag ~= 0

        % Loop by the sequence number of the concerned displacements
        for ii = 1:2

            % Calculate the dimensional-realated D_Ui_D_phi
            for jj = 1:2
                D_sK_D_phi = reshape(reshape(D_KE_D_phi_xy(:, :, jj), [], 1) * (Emin + rho_Phys(:).'.^p * (E_actual - Emin)), 64 * nelx * nely, 1);

                D_K_D_phi = sparse(iK, jK, D_sK_D_phi);
                D_K_D_phi = (D_K_D_phi + D_K_D_phi.') / 2;

                D_Ui_D_phi(jj, ii) =- adjoint_lambda(free_DOFs, ii).' * D_K_D_phi(free_DOFs, free_DOFs) * U(free_DOFs);
            end

            % Calculate the material-realated D_Ui_D_phi
            D_Ui_D_phi(3, ii) =- adjoint_lambda(free_DOFs, ii).' * D_K_D_phi_E(free_DOFs, free_DOFs) * U(free_DOFs);
            D_Ui_D_phi(4, ii) =- adjoint_lambda(free_DOFs, ii).' * D_K_D_phi_nu(free_DOFs, free_DOFs) * U(free_DOFs);

            % Calculate the loads-realated D_Ui_D_phi
            for jj = 5:4 + n_F
                temp_inv_T_sparse = sparse(F_DOFs_and_F0(:, 1), ones(size(F_DOFs_and_F0(:, 1))), inv_T{3}(:, jj - 4), 2 * (nelx + 1) * (nely + 1), 1);

                D_Ui_D_phi(jj, ii) = adjoint_lambda(free_DOFs, ii).' * temp_inv_T_sparse(free_DOFs);
            end

        end

        %% (Void) Calculate the sensitivity of compliance with respect to the standard uncertain parameters
        D_c_D_phi = zeros(4 + n_F, 1);

        % % Calculate the length-and-width-related mean compliance sensitivities
        % for jj = 1:2

        % end

        % % Calculate the Young's modulus and Poisson's ratio-related displacements sensitivities
        % for jj = 3:4

        % end

        % % Calculate the loads-related displacements sensitivities
        % for jj = 5:8

        % end

    end

    %% Calculate the sensitivities of the concerned displacements with respect to the element densities
    if nargout > 2
        node_DOF_matrix = FE_assembly.node_DOF_matrix;

        D_Ui_D_rho = zeros(nelx * nely, 2);

        for ii = 1:2
            adjoint_lambda_ii = adjoint_lambda(:, ii);

            % % The following code for calculating D_Ui_D_rho is inspired by the '88-line topology optimization' code
            D_Ui_D_rho(:, ii) = sum(-adjoint_lambda_ii(node_DOF_matrix) * KE .* U(node_DOF_matrix), 2) .* (p * (E_actual - Emin) * rho_Phys(:).^(p - 1));
        end

        %% (Void) The sensitivity of the compliance with respect to the element densities
        % ce = sum((U(node_DOF_matrix) * KE) .* U(node_DOF_matrix), 2);
        % c = sum(sum((Emin + rho_Phys.^p * (E0 - Emin)) .* ce));
        % D_c_D_rho = -p * (E0 - Emin) * rho_Phys.^(p - 1) .* ce;
        c = 0;
        D_c_D_rho = 0;
    end

end
