function [Alpha, D_Alpha_D_rho, concerned_point, inner_iter_history] = get_Alpha(...
        n_ellipsoid, Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, G_points, G_weights, ...
        rho_Phys, p, FE_assembly, BC, F_DOFs_and_F0, Ui, inv_T, fmincon_optimality_tolerance, fmincon_algorithm)
    %% This function calculates all the concerned performance values (CPVs) and the corresponding
    %  sensitivities with respect to element densities
    %
    %% ---------------------------------------------- Input ----------------------------------------------
    %   n_ellipsoid: Number of the ellipsoid convex models
    %
    %   Lx: The nominal structural length
    %   Ly: The nominal structural width
    %   nelx, nely: Element numbers in x and y directions
    %
    %   E_bounds: The bounds of the Young's modulus samples
    %   nu_bounds: The bounds of the Poisson's ratio samples
    %
    %   theta_E0, theta_nu0: The the center of the ECM of the normalized Young's modulus and Poisson's ratio uncertainties
    %
    %   G_points: Two-dimensional Gaussian points
    %   G_weights: Two-dimensional Gaussian weights (multiplied)
    %
    %   rho_Phys: The physical element density vector
    %
    %   p: The penalty factor
    %
    %   FE_assembly: The finite elements assembly information
    %
    %   BC.free_DOF: Free DOFs
    %
    %                   ┌                                ┐
    %                   |   F1_DOF    F1_nominal_value   |
    %   F_DOFs_and_F0 = │   F2_DOF    F2_nominal_value   │
    %                   │     ...           ...          │
    %                   └                                ┘
    %
    %                       The arrangement of the elements of Ui
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
    %    inv_T{3}: The inverse of the transformation matrix of the uncertain loads
    %
    %   fmincon_optimality_tolerance: The optimality tolerance of fmincon
    % 
    %   fmincon_algorithm: The algorithm adopted by fmincon
    %
    %% ---------------------------------------------- Output ----------------------------------------------
    %   Alpha: The concerned performances values
    %
    %   D_Alpha_D_rho: The sensitivities of the concerned performance values with respect to the element densities
    %
    %   concerned_point: The concerned points
    %
    %   inner_iter_history: The iteration history of the Alpha value

    %% Number of concerned displacements
    [n_Ui, ~] = size(Ui);

    %% The number of uncertain loads
    [n_F, ~] = size(F_DOFs_and_F0);

    %% Initialize the outputs
    % The concerned performance value
    Alpha = zeros(n_Ui, 1);

    % The sensitivity of the concerned performance value with respect to element densities
    D_Alpha_D_rho = zeros(nelx * nely, n_Ui);

    % The concerned point
    concerned_point = zeros(2 * n_ellipsoid, n_Ui);

    % The iteration history of Alpha in the inner loop
    inner_iter_history.Alpha = zeros(500, n_Ui);

    % The iteration steps took in the inner loop
    inner_iter_history.iter = zeros(n_Ui, 1);

    %% Calculate the Alpha values
    %  Loop by the concerned displacement responses
    for ii = 1:n_Ui
        %% Objective function
        problem.objective = @(phi) get_Alpha_objective(...
            Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, G_points, G_weights, ...
            rho_Phys, p, FE_assembly, BC, F_DOFs_and_F0, Ui, inv_T, ii, phi);

        %% Constraint function
        problem.nonlcon = @(phi) get_Alpha_constraint(phi, ii, Ui, n_ellipsoid);

        %% Initial point
        problem.x0 = zeros(4 + n_F, 1);

        %% The solver chosen
        problem.solver = 'fmincon';

        %% The options for the fmincon
        problem.options = optimoptions(...
            'fmincon', ...
            'OutputFcn', @outputFun, ...
            'OptimalityTolerance', fmincon_optimality_tolerance, ...
            'Display', 'iter', ...
            'Algorithm', fmincon_algorithm, ...
            'SpecifyObjectiveGradient', true, ...
            'SpecifyConstraintGradient', true);

        %% Set parameters for the inner loop optimization
        iteration_history.iteration_steps = 0;
        iteration_history.objective_value = [];

        %% Solve the inner loop optimization problem using 'fmincon' (Other void parameters are omitted by default.)
        [concerned_point(:, ii), Alpha(ii, 1)] = fmincon(problem);

        %% Calculate the sensitivities of the CPV with respect to the element densities
        [~, ~, ~, D_Alpha_D_rho_all] = FEA(...
            Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, G_points, G_weights, ...
            rho_Phys, p, FE_assembly, BC, F_DOFs_and_F0, Ui, inv_T, concerned_point(:, ii), 0);

        D_Alpha_D_rho(:, ii) = D_Alpha_D_rho_all(:, ii);

        %% Return the iteration history of the Alpha value
        inner_iter_history.Alpha(1:iteration_history.iter + 1, ii) = iteration_history.objective_value(:);
        inner_iter_history.iter(ii) = iteration_history.iter;

    end

    %% ------------------------------- Nested output control function -------------------------------
    function stop = outputFun(~, optimValues, state)
        stop = false;

        switch state
            case 'init'
                % hold on

            case 'iter'
                iteration_history.objective_value = [iteration_history.objective_value; optimValues.fval];
                iteration_history.iter = optimValues.iteration;

            case 'done'
                % hold off

            otherwise
        end

    end

end
