function plot_3D_to_2D_projection(W, x0, x, color, flag_isnormal)
    %% This function plots the ellipse projections of the ellipsoid

    % % Code written by Song Bai, Daming Li and Zhan Kang

    %   Input                 Description
    %  ------- ------------------------------------------
    %   W       The matrix W of the ellipsoid
    %   x0      The center of the ellipsoid
    %   x       The samples
    %   color   The color of the ellipse and the samples

    %%
    n_samples = size(x, 2);

    index_target = nchoosek(1:3, 2);

    angle_steps = linspace(0, 2 * pi, 200);

    %%
    for ii = 1:3

        % % Get the matrix W of the projection ellipse (a Schur's complement of the ECM's W)
        w_projection = get_w_projection(W, index_target(ii, :));

        x0_projection = x0(index_target(ii, :));

        % % Get the matrix A and vector b
        [Q, SIGMA] = eig(w_projection);
        A = Q * sqrt(SIGMA) * Q.';
        b = A * x0_projection;

        % % Plot the ellipse
        points_on_ellipse = A \ [cos(angle_steps) + b(1); sin(angle_steps) + b(2)];

        points_on_ellipse_projection = zeros(3, 200);
        points_on_ellipse_projection(index_target(ii, :), :) = points_on_ellipse;

        samples_projection = zeros(3, n_samples);
        samples_projection(index_target(ii, :), :) = x(index_target(ii, :), :);

        if size(color, 2) > 1 && flag_isnormal
            p = plot3(...
                points_on_ellipse_projection(1, :), ...
                points_on_ellipse_projection(2, :), ...
                points_on_ellipse_projection(3, :));

            p.LineWidth = 2;
            p.Color = color(ii);

            % % Plot the samples
            scatter3(...
                samples_projection(1, :), ...
                samples_projection(2, :), ...
                samples_projection(3, :), ...
                55, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', color(ii));

        else

            % % Plot the samples
            scatter3(...
                samples_projection(1, :), ...
                samples_projection(2, :), ...
                samples_projection(3, :), ...
                125, 'x', 'MarkerEdgeColor', color(ii), 'LineWidth', 3);
        end

    end

end
