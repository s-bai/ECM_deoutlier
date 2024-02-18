function det_lof_2d()
    %% This function implements the deterministic topology optimization that minimizing the volume fraction under displacement constraints
    %% ---------------------------------------------- Input ----------------------------------------------
    %   The JSON file for the deterministic topology optimization.
    %
    %% ---------------------------------------------- Output ----------------------------------------------
    %   Void

    %% Get execution date and time (yyyy mm dd hh mm, eg: clock_string = '202103041711')
    clock_string = get_clock();

    % Create the MAT file for saving results
    result_filename = strcat('deterministic_cantilever_beam_', clock_string, '.mat');

    %% Load the problem parameters from the JASON file
    JSON_data = jsondecode(fileread('deterministic_cantilever_beam.json'));

    % nelx, nely: Numbers of elements in x and y directions.
    nelx = JSON_data.nelx;
    nely = JSON_data.nely;

    % V0: Volume fraction of the initial design.
    V0 = JSON_data.V0;

    % p: Penalty factor.
    p = JSON_data.p;

    % rmin: Sensitivity filter radius.
    rmin = JSON_data.rmin;

    % Lx and Ly: The nominal length and width
    Lx = JSON_data.Lx;
    Ly = JSON_data.Ly;

    % E_bounds and nu_bounds: The upper and lower bounds of the Young's modulus and Poisson's ratio data samples
    E_bounds = JSON_data.E_bounds;
    nu_bounds = JSON_data.nu_bounds;

    % theta_E0, theta_nu0: The the center of the ECM of the normalized Young's modulus and Poisson's ratio samples
    theta_E0 = JSON_data.theta_E0;
    theta_nu0 = JSON_data.theta_nu0;

    %% Get the DOFs and the nominal values of the uncertain loads (N)
    [n_F, ~] = size(JSON_data.F_DOFs_and_F0);

    F_DOFs_and_F0 = zeros(n_F, 2);

    for ii = 1:n_F

        for jj = 1:2
            F_DOFs_and_F0(ii, jj) = eval(JSON_data.F_DOFs_and_F0{ii}{jj});
        end

    end

    %% Get the characteristic matrices of all uncertainties
    %    w{1}: The characteristic matrix of the uncertain length and width
    %    w{2}: The characteristic matrix of the uncertain Young's modulus and the Poisson's ratio
    %    w{3, 4, ...}: The characteristic matrices of the uncertain loads
    w = cell(1, 3);

    for ii = 1:3
        temp_variable_name = strcat('JSON_data.w', num2str(ii));
        w{ii} = eval(temp_variable_name);
    end

    %% Calculate the inverse of the transformation matrices
    %    inv_T{1}: The inverse of the transformation matrix of the uncertain length and width
    %    inv_T{2}: The inverse of the transformation matrix of the uncertain Young's modulus and the Poisson's ratio
    %    inv_T{3, 4, ...}: The inverses of the transformation matrices of the uncertain loads
    [~, inv_T] = get_inverse_T(w);

    %% Get the displacements' constraints and the corresponding reliability constraints
    %
    %            The arrangement of the elements of Ui
    %     ┌─────────┬────────────────┬──────────────────────┐
    %     │    Ui   │ Ui_upper_limit │   eta_lower_limit    │
    %     ├─────────┼────────────────┼──────────────────────┤
    %     │   DOF1  │     U_DOF1     │   eta_lower_limit1   │
    %     │   DOF2  │     U_DOF2     │   eta_lower_limit2   │
    %     │   ...   │     ...        │          ...         │
    %     └─────────┴────────────────┴──────────────────────┘
    %
    %     E.g.
    %     ┌──────────────────────────────┬────────────────┬──────────────────────┐
    %     │            Ui                │ Ui_upper_limit │   eta_lower_limit    │
    %     ├──────────────────────────────┼────────────────┼──────────────────────┤
    %     │   (nelx + 2) * (nely + 1)    │      -1.0      │          1.0         │
    %     │ 2 * (nelx + 1) * (nely + 1)  │      -2.2      │          1.0         │
    %     └──────────────────────────────┴────────────────┴──────────────────────┘
    Ui = zeros(n_F, 3);

    for ii = 1:2

        for jj = 1:3
            Ui(ii, jj) = eval(JSON_data.Ui{ii}{jj});
        end

    end

    %% Calculate two-dimensional Gaussian points matrix
    [G_points, G_weights] = get_G_points(3);

    %% Numbers of constraints and design variables
    [n_constraint, ~] = size(Ui);
    n_design_variable = nelx * nely;

    %% Boundary conditions of the cantilever beam with two loads
    % Fixed DOF
    BC.fixed_DOFs = eval(JSON_data.fixed_DOFs{1});

    % Free DOF
    BC.free_DOFs = setdiff(1:2 * (nelx + 1) * (nely + 1), BC.fixed_DOFs);

    %% Parameters for the MMA optimizer
    [rho_min, rho_max, ...
            rho_old1_vector, rho_old2_vector, ...
            lower_asymptotes, upper_asymptotes, ...
            a0_MMA, a_MMA, c_MMA, d_MMA] = get_MMA_parameter(nelx, nely, n_constraint);

    constraint_displacements = zeros(2, 1);

    %% Calculate the filter parameters
    [H, Hs] = get_filter(nelx, nely, rmin);

    %% Calculate finite element nodes information
    FE_assembly = get_FE_assembly(nelx, nely);

    %% Iteration step counter and convergence criteria
    iter = 0;
    change = 1;
    max_change = JSON_data.max_change;
    max_iter_step = JSON_data.max_iter_step;

    %% Initialize design variables
    rho = repmat(V0, nely, nelx);
    rho_Phys = rho;

    rho_Phys_history = zeros(nely, nelx, max_iter_step + 1);
    rho_Phys_history(:, :, iter + 1) = rho_Phys;

    iter_histories = zeros(max_iter_step + 1, 5);

    % Set all the uncertain parameters to be their nominal values
    phi = zeros(6, 1);

    %% Scaling factor of the objective function
    % (Used to keep the objective function and the constraint functions' values at the same order of magnitude)
    objective_epsilon = JSON_data.objective_epsilon;

    %% Open figure window
    figure('Name', 'Deterministic design of the cantilever beam', 'NumberTitle', 'off');

    % Start the stopwatch timer
    t_start = tic;

    %% Deterministic topology optimization
    while (change >= max_change) && (iter <= max_iter_step)
        iter = iter + 1;

        %% Calculate the objective function (the total volume fraction)
        objective_volume = objective_epsilon * mean(mean(rho));
        % objective_volume = objective_epsilon * sum(sum(rho));

        % The 1st and 2nd order sensitivities of the objective function
        D_objective_volume = objective_epsilon * ones(nelx * nely, 1) / (nelx * nely);
        % D_objective_volume = objective_epsilon * ones(nelx * nely, 1);
        D2_objective_volume = 0 * D_objective_volume;

        %% Calculate the constraint functions (the concerned displacements)
        [U, ~, ~, DUi_D_rho] = FEA(...
            Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, G_points, G_weights, ...
            rho_Phys, p, FE_assembly, BC, F_DOFs_and_F0, Ui, inv_T, phi, 0);

        constraint_displacements(:, 1) = -U(Ui(:, 1)) + Ui(:, 2);

        % The 1st and 2nd order sensitivities of the constraint functions
        D_constraint_displacements = -DUi_D_rho;

        %% Filter the sensitivities of the objective function and the constraint functions
        % The density filter
        D_objective_volume = H * (D_objective_volume ./ Hs);

        for ii = 1:2
            D_constraint_displacements(:, ii) = (H * (D_constraint_displacements(:, ii) ./ Hs));
        end

        D_constraint_displacements = D_constraint_displacements.';

        D2_constraint_displacements = 0 * D_constraint_displacements;

        %% Update the design variables using the MMA optimizer
        rho_old = rho(:);

        [rho_new_vector, ~, ~, ~, ~, ~, ~, ~, ~, lower_asymptotes, upper_asymptotes] = mmasub(...
            n_constraint, n_design_variable, iter, ...
            rho(:), rho_min(:), rho_max(:), rho_old1_vector, rho_old2_vector, ...
            objective_volume, D_objective_volume(:), D2_objective_volume(:), ...
            constraint_displacements, D_constraint_displacements, D2_constraint_displacements, ...
            lower_asymptotes, upper_asymptotes, a0_MMA, a_MMA, c_MMA, d_MMA);

        % % Filter the densities by the density filter
        rho_Phys = (H * rho_new_vector) ./ Hs;
        rho_Phys = reshape(rho_Phys, nely, nelx);

        rho = reshape(rho_new_vector, nely, nelx);

        rho_old2_vector = rho_old1_vector;
        rho_old1_vector = rho_old;

        %% Calculate the change between the two consecutive iteration steps
        change = max(abs(rho_new_vector - rho_old));

        %% Print optimization results to the command window
        fprintf('Iter. #: %d     V_fraction: %.2f     U1: %.2f     U2: %.2f    change: %.3f\n', ...
            iter, objective_volume / objective_epsilon, U(Ui(1, 1)), U(Ui(2, 1)), change);

        %% Plot densities
        colormap(gray);
        imagesc(1 - rho_Phys);
        caxis([0 1]);
        axis equal;
        axis off;
        drawnow;

        %% Save iteration histories to the result file
        rho_Phys_history(:, :, iter + 1) = rho_Phys;
        iter_histories(iter, :) = [iter, objective_volume / objective_epsilon, U(Ui(1, 1)), U(Ui(2, 1)), change];

        save(result_filename, 'iter', 'rho_Phys_history', 'iter_histories');

    end

    % Stop the stopwatch timer
    total_runtime = toc(t_start);

    save(result_filename, 'total_runtime', '-append');

end
