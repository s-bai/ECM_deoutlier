function MCS_u()
    %% This function implements the Monte Carlo simulation to verify the concerned displacements' feasiability under the constraints
    % ----------------------------------------------- Input -----------------------------------------------
    %
    % ----------------------------------------------- Output -----------------------------------------------
    %

    %% Load the problem parameters from the JASON file
    JSON_data = jsondecode(fileread('cantilever_beam_2D.json'));

    % nelx, nely: Numbers of elements in x and y directions.
    nelx = JSON_data.nelx;
    nely = JSON_data.nely;

    % V0: Volume fraction of the initial design.
    result_filename = 'deterministic_cantilever_beam.mat';
    result_data = load(result_filename);

    rho_Phys = result_data.rho_Phys_history(:, :, result_data.iter + 1);

    %% Create result file for the Monte Carlo simulation results
    MCS_result_filename = 'MCS_cantilever_beam_2D_.mat';

    % p: Penalty factor.
    p = JSON_data.p;

    % Lx and Ly: The nominal length and width (in mm)
    Lx = JSON_data.Lx;
    Ly = JSON_data.Ly;

    % E_bounds and nu_bounds: The upper and lower bounds of the Young's modulus (in N/mm^2) and Poisson's ratio data samples
    E_bounds = JSON_data.E_bounds;
    nu_bounds = JSON_data.nu_bounds;

    % theta_E0, theta_nu0: The the center of the ECM of the normalized Young's modulus and Poisson's ratio uncertainties
    theta_E0 = JSON_data.theta_E0;
    theta_nu0 = JSON_data.theta_nu0;

    %% Get the DOFs and the nominal values of the uncertain loads (N)
    F_DOFs_and_F0 = zeros(4, 2);

    for ii = 1:4

        for jj = 1:2
            F_DOFs_and_F0(ii, jj) = eval(JSON_data.F_DOFs_and_F0{ii}{jj});
        end

    end

    %% Get the characteristic matrices of all uncertainties
    %    w{1}: The characteristic matrix of the uncertain length and width
    %    w{2}: The characteristic matrix of the uncertain Young's modulus and the Poisson's ratio
    %    w{3, 4, ...}: The characteristic matrices of the uncertain loads
    w = cell(1, 4);

    for ii = 1:4
        temp_variable_name = strcat('JSON_data.w', num2str(ii));
        w{ii} = eval(temp_variable_name);
    end

    %% Calculate the inverse of the transformation matrices
    %    inv_T{1}: The inverse of the transformation matrix of the uncertain length and width
    %    inv_T{2}: The inverse of the transformation matrix of the uncertain Young's modulus and the Poisson's ratio
    %    inv_T{3, 4, ...}: The inverses of the transformation matrices of the uncertain loads
    [n_ellipsoid, inv_T] = get_inverse_T(w);

    %% Calculate the parameters for the uncertain loads
    F_parameter = get_F_parameter(nelx, nely, inv_T, F_DOFs_and_F0);

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
    Ui = zeros(2, 3);

    for ii = 1:2

        for jj = 1:3
            Ui(ii, jj) = eval(JSON_data.Ui{ii}{jj});
        end

    end

    %% Get the samples number for the Monte Carlo simulation
    n_MCS_samples = JSON_data.n_MCS_samples;
    Ui_MCS_results = zeros(2, n_MCS_samples);

    %% Calculate two-dimensional Gaussian points matrix
    [G_points, G_weights] = get_G_points(3);

    %% --------------------- Boundary conditions of the cantilever beam with two loads ---------------------
    % All degrees of freedom
    BC.all_DOFs = 1:2 * (nelx + 1) * (nely + 1);

    % Fixed DOF
    BC.fixed_DOFs = 1:2 * (nely + 1);

    % Free DOF
    BC.free_DOFs = setdiff(BC.all_DOFs, BC.fixed_DOFs);

    %% Calculate finite element nodes information
    FE_assembly = get_FE_assembly(nelx, nely);

    %% Generate Monte Carlo simulation samples
    phi_MCS_samples = get_phi_MCS_samples(n_ellipsoid, n_MCS_samples, 0);

    %% Implement the Monte Carlo simulation
    parfor ii = 1:n_MCS_samples
        [U] = FEA(Lx, Ly, nelx, nely, E_bounds, nu_bounds, theta_E0, theta_nu0, G_points, G_weights, rho_Phys, p, FE_assembly, BC, F_parameter, Ui, inv_T, phi_MCS_samples(:, ii), 0);
        Ui_MCS_results(:, ii) = U(Ui(:, 1));
    end

    %%
    fprintf('The maximum concerned displacements are %.2f and %.2f\n', min(Ui_MCS_results(1, :)), min(Ui_MCS_results(2, :)));

    %% Save the Monte Carlo simulation results to the file
    save(MCS_result_filename, 'Ui_MCS_results');

end
