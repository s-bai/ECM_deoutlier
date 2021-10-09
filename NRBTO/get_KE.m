function [KE, D_KE_D_phi_xy, D_KE_D_E, D_KE_D_nu] = ...
        get_KE(Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, phi, inv_T, G_points, G_weights)
    %% This function calculates the element stiffness matrix KE and its sensitivities with respect to the normalized dimensional and material uncertainties
    %
    %  Note that the element stiffness matrix is not multiplied by the Young's modulus
    %
    %% ---------------------------------- Input ----------------------------------
    %   Lx: The nominal structural length
    %   Ly: The nominal structural width
    %   nelx, nely: Element numbers in the x and y directions
    %
    %   E_bounds: The bounds of the Young's modulus samples
    %   nu_bounds: The bounds of the Poisson's ratio samples
    %
    %   theta_E0, theta_nu0: The the center of the ECM of the normalized Young's modulus and Poisson's ratio uncertainties
    %
    %   phi: Normalized uncertain parameters vector
    %       phi(1), phi(2): For the length and width uncertainties
    %       phi(3), phi(4): For the Young's modulus and the Poisson's ratio uncertainties
    %       phi(5), phi(6): For the loads uncertainties of F1
    %       phi(7), phi(8): For the loads uncertainties of F2
    %
    %   inv_T: The inverse of the transformation matrix
    %    inv_T{1}: The inverse of the transformation matrix of the uncertain length and width
    %    inv_T{2}: The inverse of the transformation matrix of the uncertain Young's modulus and the Poisson's ratio
    %    inv_T{3, 4, ...}: The inverses of the transformation matrices of the uncertain loads
    %
    %   G_points: Two-dimensional Gaussian points
    %   G_weights: Two-dimensional Gaussian weights (The tensor product of the weights)
    %
    %% ---------------------------------- Output ----------------------------------
    %   KE: The element stiffness matrix
    %   D_KE_D_phi_xy(:, :, 1): The sensitivities of KE with respect to the normalized length uncertainty
    %   D_KE_D_phi_xy(:, :, 2): The sensitivities of KE with respect to the normalized width uncertainty
    %
    %   D_KE_D_E: The sensitivities of KE with respect to the Young's modulus
    %   D_KE_D_nu: The sensitivities of KE with respect to the Poisson's ratio (without multiplying the E)
    %
    %   Note that the sensitivities of KE with respect to the normalized parameters of uncertain Young's modulus and Poisson's ratio are calculated in the FEA subroutine

    %% Initialize the 3D matrix storing the matrix B and the sensitivities of B at all quadrature points
    G_points_number = max(size(G_points));
    matrix_B_quadrature = zeros(3, 8, G_points_number);

    if nargout > 1
        Dmatrix_B_quadrature_D_phi_x = zeros(3, 8, G_points_number);
        Dmatrix_B_quadrature_D_phi_y = zeros(3, 8, G_points_number);

        [matrix_D, Dmatrix_D_D_E, Dmatrix_D_D_nu] = get_matrix_D(E_bounds, nu_bounds, theta_E0, theta_nu0, phi, inv_T);

        D_KE_D_phi_xy = zeros(8, 8, 2);

    else
        [matrix_D] = get_matrix_D(E_bounds, nu_bounds, theta_E0, theta_nu0, phi, inv_T);

    end

    %% Calculate the actual element length and width
    [lxx, lyy] = get_lxx_lyy(Lx, Ly, nelx, nely, phi, inv_T);

    %% Calculate the Jacobian determinant
    Jacobian_determinant = lxx * lyy / 4;

    if nargout > 1
        D_Jacobian_determinant_D_phi_x = 1/4 * (inv_T{1}(1, 1) / nelx * lyy + lxx * inv_T{1}(2, 1) / nely);
        D_Jacobian_determinant_D_phi_y = 1/4 * (inv_T{1}(1, 2) / nelx * lyy + lxx * inv_T{1}(2, 2) / nely);
    end

    %% Calculate the values of matrix B and matrix D at all quadrature points
    for ii = 1:G_points_number

        if nargout > 1
            [matrix_B_quadrature(:, :, ii), Dmatrix_B_quadrature_D_phi_x(:, :, ii), Dmatrix_B_quadrature_D_phi_y(:, :, ii)] = ...
                get_matrix_B(G_points(ii, :), lxx, lyy, nelx, nely, inv_T);

        else
            [matrix_B_quadrature(:, :, ii)] = ...
                get_matrix_B(G_points(ii, :), lxx, lyy, nelx, nely, inv_T);
        end

    end

    %% Calculate KE and the sensitivities of KE with respect to phi_x, phi_y, E and nu
    KE = zeros(8, 8);

    if nargout > 1
        D_KE_D_phi_x = zeros(8, 8);
        D_KE_D_phi_y = zeros(8, 8);

        D_KE_D_E = zeros(8, 8);
        D_KE_D_nu = zeros(8, 8);

    end

    for ii = 1:G_points_number
        KE = KE + ...
            matrix_B_quadrature(:, :, ii).' * matrix_D * matrix_B_quadrature(:, :, ii) * Jacobian_determinant * G_weights(ii);

        if nargout > 1
            D_KE_D_phi_x = D_KE_D_phi_x + ...
                Dmatrix_B_quadrature_D_phi_x(:, :, ii).' * matrix_D * matrix_B_quadrature(:, :, ii) * Jacobian_determinant * G_weights(ii) + ...
                matrix_B_quadrature(:, :, ii).' * matrix_D * Dmatrix_B_quadrature_D_phi_x(:, :, ii) * Jacobian_determinant * G_weights(ii) + ...
                matrix_B_quadrature(:, :, ii).' * matrix_D * matrix_B_quadrature(:, :, ii) * D_Jacobian_determinant_D_phi_x * G_weights(ii);

            D_KE_D_phi_y = D_KE_D_phi_y + ...
                Dmatrix_B_quadrature_D_phi_y(:, :, ii).' * matrix_D * matrix_B_quadrature(:, :, ii) * Jacobian_determinant * G_weights(ii) + ...
                matrix_B_quadrature(:, :, ii).' * matrix_D * Dmatrix_B_quadrature_D_phi_y(:, :, ii) * Jacobian_determinant * G_weights(ii) + ...
                matrix_B_quadrature(:, :, ii).' * matrix_D * matrix_B_quadrature(:, :, ii) * D_Jacobian_determinant_D_phi_y * G_weights(ii);

            D_KE_D_E = D_KE_D_E + ...
                matrix_B_quadrature(:, :, ii).' * Dmatrix_D_D_E * matrix_B_quadrature(:, :, ii) * Jacobian_determinant * G_weights(ii);

            D_KE_D_nu = D_KE_D_nu + ...
                matrix_B_quadrature(:, :, ii).' * Dmatrix_D_D_nu * matrix_B_quadrature(:, :, ii) * Jacobian_determinant * G_weights(ii);

        end

    end

    if nargout > 1
        D_KE_D_phi_xy(:, :, 1) = D_KE_D_phi_x;
        D_KE_D_phi_xy(:, :, 2) = D_KE_D_phi_y;
    end

end

% ------------------------------------- Local functions -------------------------------------

function [matrix_B, Dmatrix_B_D_phi_x, Dmatrix_B_D_phi_y] = get_matrix_B(xi_eta, lxx, lyy, nelx, nely, inv_T)
    %% This local function calculates matrix B
    %
    %% Input:
    %   xi_eta = [xi, eta]
    %   lxx, lyy: The actual element length and width
    %   nelx, nely: The element numbers in the length and width directions
    %
    %   inv_T: See the header comments
    %
    %% Output:
    %   matrix_B: The matrix B
    %   Dmatrix_B_D_phi_x: Sensitivity of matrix B with respect to the normalized length uncertainty
    %   Dmatrix_B_D_phi_y: Sensitivity of matrix B with respect to the normalized width uncertainty

    %%
    xi = xi_eta(1);
    eta = xi_eta(2);

    %% Matrix B
    matrix_B = [-(1 - eta) / (2 * lxx), 0, (1 - eta) / (2 * lxx), 0, (1 + eta) / (2 * lxx), 0, -(1 + eta) / (2 * lxx), 0; ...
                0, -(1 - xi) / (2 * lyy), 0, -(1 + xi) / (2 * lyy), 0, (1 + xi) / (2 * lyy), 0, (1 - xi) / (2 * lyy); ...
                -(1 - xi) / (2 * lyy), -(1 - eta) / (2 * lxx), -(1 + xi) / (2 * lyy), (1 - eta) / (2 * lxx), (1 + xi) / (2 * lyy), (1 + eta) / (2 * lxx), (1 - xi) / (2 * lyy), -(1 + eta) / (2 * lxx)];

    %% The sensitivity of matrix B with respect to phi_x and phi_y
    if nargout > 1
        Dmatrix_B_D_lxx = zeros(3, 8);
        Dmatrix_B_D_lyy = zeros(3, 8);

        Dmatrix_B_D_lxx(1, :) = [(1 - eta) / (2 * lxx^2), 0, -(1 - eta) / (2 * lxx^2), 0, -(1 + eta) / (2 * lxx^2), 0, (1 + eta) / (2 * lxx^2), 0];

        Dmatrix_B_D_lxx(3, :) = [0, (1 - eta) / (2 * lxx^2), 0, -(1 - eta) / (2 * lxx^2), 0, -(1 + eta) / (2 * lxx^2), 0, (1 + eta) / (2 * lxx^2)];

        Dmatrix_B_D_lyy(2, :) = [0, (1 - xi) / (2 * lyy^2), 0, (1 + xi) / (2 * lyy^2), 0, -(1 + xi) / (2 * lyy^2), 0, -(1 - xi) / (2 * lyy^2)];

        Dmatrix_B_D_lyy(3, :) = [(1 - xi) / (2 * lyy^2), 0, (1 + xi) / (2 * lyy^2), 0, -(1 + xi) / (2 * lyy^2), 0, -(1 - xi) / (2 * lyy^2), 0];

        Dmatrix_B_D_phi_x = Dmatrix_B_D_lxx * inv_T{1}(1, 1) / nelx + Dmatrix_B_D_lyy * inv_T{1}(2, 1) / nely;
        Dmatrix_B_D_phi_y = Dmatrix_B_D_lxx * inv_T{1}(1, 2) / nelx + Dmatrix_B_D_lyy * inv_T{1}(2, 2) / nely;
    end

end

function [matrix_D, Dmatrix_D_D_E, Dmatrix_B_D_nu] = get_matrix_D(E_bounds, nu_bounds, theta_E0, theta_nu0, phi, inv_T)
    %% This local function calculates the matrix D
    % Input:
    %   See the header comments 
    %
    %% Output:
    %   matrix_D: The matrix D
    %   Dmatrix_D_D_E: The sensitivity of matrix D with respect to the uncertain Young's modulus
    %   Dmatrix_B_D_nu: The sensitivity of matrix D with respect to the uncertain Poisson's ratio

    %% Calculate the actual Poisson's ratio
    [~, nu_actual] = get_E_nu(E_bounds, nu_bounds, theta_E0, theta_nu0, phi, inv_T);

    %% Calculate the matrix D (without multiplying the E)
    matrix_D = ...
        1 / (1 - nu_actual^2) * ...
        [1, nu_actual, 0; ...
            nu_actual, 1, 0; ...
            0, 0, (1 - nu_actual) / 2];

    %% Calculate the sensitivity of matrix D with respect to E
    Dmatrix_D_D_E = ...
        1 / (1 - nu_actual^2) * ...
        [1, nu_actual, 0; ...
            nu_actual, 1, 0; ...
            0, 0, (1 - nu_actual) / 2];

    %% Calculate the sensitivity of matrix D with respect to nu (without multiplying the E)
    if nargout > 1

        Dmatrix_B_D_nu = zeros(3, 3);

        Dmatrix_B_D_nu(1) = 2 * nu_actual;
        Dmatrix_B_D_nu(2) = 1 + nu_actual^2;
        Dmatrix_B_D_nu(5) = 2 * nu_actual;
        Dmatrix_B_D_nu(9) = -1/2 * (1 - nu_actual)^2;
        Dmatrix_B_D_nu = Dmatrix_B_D_nu.' + Dmatrix_B_D_nu - diag(diag(Dmatrix_B_D_nu));
        Dmatrix_B_D_nu = 1 / ((1 - nu_actual^2)^2) * Dmatrix_B_D_nu;

    end

end

function [lxx, lyy] = get_lxx_lyy(Lx, Ly, nelx, nely, phi, inv_T)
    %% This local function calculates the actual structure length and width
    %
    %% Input:
    %   See the header comments
    %
    %% Output:
    %   lxx, lyy: The actual element length and width

    %%
    lxx_lyy = inv_T{1} * [phi(1); phi(2)] + [Lx; Ly];

    lxx = lxx_lyy(1) / nelx;
    lyy = lxx_lyy(2) / nely;

end
