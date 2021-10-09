function nrbto_lof_2d(JSON_filename)
    %% ----------- Non-probabilistic Reliability-Based Topology Optimization of 2D structure -----------
    %
    %
    %                              Song Bai, Daming Li and Zhan Kang
    %
    %      Dept. of Engineering Mechanics, Dalian University of Technology, Dalian, China, 116024
    %
    %
    %
    %% ---------------------------------------------- Input ----------------------------------------------
    %   The JSON file.
    %
    %% ---------------------------------------------- Output ---------------------------------------------
    %   Results will be printed to the screen and saved to the file.

    %% -------------------------------- Load parameters from the JSON file -------------------------------
    if nargin > 0
        JSON_data = jsondecode(fileread(JSON_filename));
    else
        JSON_data = jsondecode(fileread('NRBTO_cantilever_beam.json'));
    end

    % nelx, nely: Numbers of elements in x and y directions.
    nelx = JSON_data.nelx;
    nely = JSON_data.nely;

    % V0: Volume fraction of the initial design.
    V0 = JSON_data.V0;

    % p: The penalty factor.
    p = JSON_data.p;

    % rmin: Sensitivity filter radius.
    rmin = JSON_data.rmin;

    % Lx and Ly: The nominal length and width
    Lx = JSON_data.Lx;
    Ly = JSON_data.Ly;

    % E_bounds and nu_bounds: The upper and lower bounds of the Young's modulus and Poisson's ratio data samples
    E_bounds = JSON_data.E_bounds;
    nu_bounds = JSON_data.nu_bounds;

    % theta_E0, theta_nu0: The the center of the ECM of the normalized Young's modulus and Poisson's ratio uncertainties
    theta_E0 = JSON_data.theta_E0;
    theta_nu0 = JSON_data.theta_nu0;

    %% Get the DOFs and the nominal values of the uncertain loads
    % Get the number of uncertain loads, which is also the number of concerned displacements
    [n_F, ~] = size(JSON_data.F_DOFs_and_F0);
    n_Ui = n_F;

    F_DOFs_and_F0 = zeros(n_F, 2);

    for ii = 1:n_F

        for jj = 1:2
            F_DOFs_and_F0(ii, jj) = eval(JSON_data.F_DOFs_and_F0{ii}{jj});
        end

    end

    %% Get the normalized characteristic matrices of all uncertainties
    %    w{1}: The characteristic matrix of the uncertain length and width
    %    w{2}: The characteristic matrix of the uncertain Young's modulus and the Poisson's ratio
    %    w{3}: The characteristic matrices of all uncertain loads
    w = cell(1, 3);

    for ii = 1:3
        temp_variable_name = strcat('JSON_data.w', num2str(ii));
        w{ii} = eval(temp_variable_name);
    end

    %% Calculate the inverse of the transformation matrices
    %    inv_T{1}: The inverse of the transformation matrix of the uncertain length and width
    %    inv_T{2}: The inverse of the transformation matrix of the uncertain Young's modulus and the Poisson's ratio
    %    inv_T{3}: The inverse of the transformation matrix of the uncertain loads
    [n_ellipsoid, inv_T] = get_inverse_T(w);

    %% Get the displacements' constraints and the corresponding reliability constraints (See the README for details)
    Ui = zeros(n_F, 3);

    for ii = 1:2

        for jj = 1:3
            Ui(ii, jj) = eval(JSON_data.Ui{ii}{jj});
        end

    end

    %% Calculate two-dimensional Gaussian quadrature points matrix
    [G_points, G_weights] = get_G_points(3);

    %% Get the numbers of constraints and design variables
    [n_constraint, ~] = size(Ui);
    n_design_variable = nelx * nely;

    %% Boundary conditions
    % Fixed DOFs
    BC.fixed_DOFs = eval(JSON_data.fixed_DOFs{1});

    % Free DOFs
    BC.free_DOFs = setdiff(1:2 * (nelx + 1) * (nely + 1), BC.fixed_DOFs);

    %% --------------------- Parameters of the MMA optimizer ---------------------
    [rho_min, rho_max, rho_old1_vector, rho_old2_vector, ...
            lower_asymptotes, upper_asymptotes, a0_MMA, a_MMA, c_MMA_base, d_MMA] = get_MMA_parameter(nelx, nely, n_constraint);

    % Rescale the magnitude of c_MMA
    c_MMA = JSON_data.c_MMA_scaling * c_MMA_base;

    %% Get the filter parameters
    [H, Hs] = get_filter(nelx, nely, rmin);

    %% Calculate finite element nodes information
    FE_assembly = get_FE_assembly(nelx, nely);

    %% Initialize design variables
    rho = repmat(V0, nely, nelx);
    rho_Phys = rho;

    %% Initialize the iteration step counter and convergence criteria
    iter = 0;
    change = 1;
    max_change = JSON_data.max_change;
    max_iter_step = JSON_data.max_iter_step;

    %% Initilize results storage arrays
    rho_history = zeros(nely, nelx, max_iter_step + 1);
    rho_Phys_history = zeros(nely, nelx, max_iter_step + 1);
    volume_fraction_history = zeros(1, max_iter_step + 1);
    concerned_point_history = zeros(2 * n_ellipsoid, n_constraint, max_iter_step + 1);
    Alpha_history = zeros(n_Ui, max_iter_step + 1);
    inner_Alpha_history = zeros(500, n_Ui, max_iter_step + 1);
    inner_iter_steps = zeros(n_Ui, max_iter_step + 1);
    y_limits = zeros(1, 2);
    plot_handle_concerned_point = cell(3, n_Ui);

    figure_number_top_evo = 0;
    figure_number_iter_history = 0;
    figure_number_CPts = zeros(1, n_Ui);

    %% Scaling factor of the objective function
    %   (Used to keep the objective function and the constraint functions' values at the same order of magnitude)
    objective_epsilon = JSON_data.objective_epsilon;

    %% The optimality tolerance of fmincon
    fmincon_optimality_tolerance = JSON_data.fmincon_optimality_tolerance;

    %% --------------------------------- The start of the optimization process ---------------------------------

    % Get the execution date and time (yyyy mm dd hh mm, eg: clock_string = '202103041711')
    clock_string = get_clock();

    % Create the MAT file for saving results
    result_file = strcat('NRBTO_', clock_string, '.mat');

    % Start the stopwatch timer
    t_start = tic;

    while (change >= max_change) && (iter <= max_iter_step)
        iter = iter + 1;

        %% Increase penalty factor by 0.2 every 10 steps after the 100th step
        if JSON_data.flag_increase_P

            if iter > 100 && mod(iter, 10) == 0
                p = p + 0.2;
            end

        end

        %% Calculate the NRBTO objective function (the total volume fraction)
        volume_fraction_history(1, iter) = mean(mean(rho));
        objective_volume = objective_epsilon * mean(mean(rho));

        % % The 1st and 2nd order sensitivities of the objective function with respect to the design variables
        D_objective_volume = objective_epsilon * ones(nelx * nely, 1) / (nelx * nely);
        D2_objective_volume = 0 * D_objective_volume;

        %% Calculate the NRBTO constraint functions
        %       Alpha: n_Ui-by-1 vector
        %       D_Alpha_D_rho: (nelx * nely)-by-n_Ui matrix
        %       concerned_point: n_phi-by-n_Ui matrix
        %       inner_iter_history.Alpha: 500-by-n_Ui matrix
        %       inner_iter_steps: n_Ui-by-1 vector
        [Alpha, D_Alpha_D_rho, concerned_point, inner_iter_history] = get_Alpha(...
            n_ellipsoid, Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, G_points, G_weights, ...
            rho_Phys, p, FE_assembly, BC, F_DOFs_and_F0, Ui, inv_T, fmincon_optimality_tolerance, JSON_data.fmincon_algorithm{1});

        % % Notice that the Alphas in the paper are larger than zero, thus one needs to add minus sign before feed them to fmincon
        constraint_Alpha = -Alpha;

        Alpha_history(:, iter) = Alpha(:, 1);

        concerned_point_history(:, :, iter) = concerned_point;

        inner_Alpha_history(:, :, iter) = inner_iter_history.Alpha(:, :);

        inner_iter_steps(:, iter) = inner_iter_history.iter(:);

        %% Apply density filter to the sensitivities
        D_objective_volume = H * (D_objective_volume(:) ./ Hs);

        D_constraint_Alpha = zeros(nelx * nely, 2);

        % % Add minus sign to the constraints sensitivities accordingly
        for ii = 1:n_Ui
            D_constraint_Alpha(:, ii) =- H * (D_Alpha_D_rho(:, ii) ./ Hs);
        end

        % Transpose the sensitivities matrix of the constraint functions
        D_constraint_Alpha = D_constraint_Alpha.';

        % The 2nd order sensitivities of the constraint functions
        D2_constraint_Alpha = 0 * D_constraint_Alpha;

        %% Update the design variables using the MMA
        rho_old = rho(:);

        [rho_new_vector, ~, ~, ~, ~, ~, ~, ~, ~, lower_asymptotes, upper_asymptotes] = mmasub(...
            n_constraint, n_design_variable, iter, ...
            rho(:), rho_min(:), rho_max(:), rho_old1_vector, rho_old2_vector, ...
            objective_volume, D_objective_volume(:), D2_objective_volume(:), ...
            constraint_Alpha, D_constraint_Alpha, D2_constraint_Alpha, ...
            lower_asymptotes, upper_asymptotes, a0_MMA, a_MMA, c_MMA, d_MMA);

        rho = reshape(rho_new_vector, nely, nelx);
        rho_history(:, :, iter) = rho;

        % % Calculate the physical element densities
        rho_Phys = reshape(((H * rho(:)) ./ Hs), nely, nelx);
        rho_Phys_history(:, :, iter) = rho_Phys;

        rho_old2_vector = rho_old1_vector;
        rho_old1_vector = rho_old;

        %% Calculate the difference between two consecutive iteration steps
        change = max(abs(rho_new_vector - rho_old));

        %% Plot the topology layout evolution to Figure #1
        if iter == 1
            % % Create figure window for the topology layout evolution (Figure #1)
            figure_number_top_evo = figure('Name', 'Topology layout evolution.');
        end

        figure(figure_number_top_evo);
        colormap(gray);
        imagesc(1 - rho_Phys);
        caxis([0 1]);
        axis equal;
        axis off;
        drawnow;

        %% Plot the iteration progress of the objective (the volume fraction) to Figure #2
        if iter == 1
            y_limits(1) = volume_fraction_history(1, iter);
            y_limits(2) = volume_fraction_history(1, iter);
        else

            if volume_fraction_history(1, iter) > y_limits(2)
                y_limits(2) = volume_fraction_history(1, iter);
            end

            if volume_fraction_history(1, iter) < y_limits(1)
                y_limits(1) = volume_fraction_history(1, iter);
            end

        end

        %%
        if iter == 1
            figure_number_iter_history = figure('Name', 'V_fraction history');
        end

        get_iteration_history_visualization(figure_number_iter_history, iter, volume_fraction_history, y_limits);

        %% Plot the concerned points to Figure #3 and after
        if iter == 1

            for ii = 1:n_Ui
                figure_number_CPts(ii) = figure('Name', strcat('The CPts of U#', num2str(ii)));
            end

        end

        plot_handle_concerned_point = ...
            get_concerned_point_visualization(figure_number_CPts, iter, plot_handle_concerned_point, concerned_point, 3, Ui);

        %% Save iteration histories to the result file
        if isfile(result_file)
            save(result_file, ...
                'volume_fraction_history', 'rho_history', 'rho_Phys_history', 'Alpha_history', ...
                'concerned_point_history', 'inner_Alpha_history', 'inner_iter_steps', '-append');
        else
            save(result_file, ...
                'volume_fraction_history', 'rho_history', 'rho_Phys_history', 'Alpha_history', ...
                'concerned_point_history', 'inner_Alpha_history', 'inner_iter_steps');
        end

        %% Print optimization results to the command window
        fprintf('It.:%5i  V_fraction:%.4f  change:%7.3f\n', iter, objective_volume / objective_epsilon, change);
        fprintf('The alpha values of the concerned displacements are:\n');

        for ii = 1:2
            fprintf('%.3f   ', Alpha_history(ii, iter));
        end

        fprintf('\n\n');
    end

    % Stop the stopwatch timer
    total_runtime = toc(t_start);
    % % ---------------------------------- The end of the optimization process ----------------------------------

    %% Print the termination status to the command line window
    if iter < max_iter_step
        fprintf('\nSolution converged after %d iteration steps.\n', iter);

    else
        fprintf('\nIteration terminated after maximum allowable iteration steps.\n');
        fprintf('The change between the last two iteration steps is %.2f\n', change);

    end

    %% Save results to the file
    save(result_file, 'iter', 'total_runtime', '-append');

end
