function phi_MCS_samples = get_phi_MCS_samples(n_ellipsoid, n_MCS_samples, plot_flag)
    %% This function generates the Monte Carlo simulation samples of the normalized uncertain parameters
    % ----------------------------------------------- Input -----------------------------------------------
    % n_ellipsoid: The number of the ellipsoid convex models
    % n_MCS_samples: The number of the Monte Carlo simulation samples
    %
    %% ----------------------------------------------- Output -----------------------------------------------
    % phi_MCS_samples: The Monte Carlo simulation samples of the normalized uncertain parameters (stored column-wisely)
    %

    %%
    phi_MCS_samples = zeros(2 * n_ellipsoid, n_MCS_samples);

    for ii = 1:n_ellipsoid
        temp_theta = 2 * pi * rand(1, n_MCS_samples);
        temp_radius = sqrt(rand(1, n_MCS_samples));

        phi_MCS_samples(2 * ii - 1, :) = temp_radius .* cos(temp_theta);
        phi_MCS_samples(2 * ii, :) = temp_radius .* sin(temp_theta);
    end

    %%
    if plot_flag ~= 0

        for ii = 1:n_ellipsoid
            figure_title = strcat('The MC samples in ellipsoid #', num2str(ii));
            figure('Name', figure_title, 'NumberTitle', 'off');

            scatter(phi_MCS_samples(2 * ii - 1, :), phi_MCS_samples(2 * ii, :), '.');
            pbaspect([1 1 1]);
        end

    end
