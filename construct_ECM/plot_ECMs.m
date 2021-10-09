function plot_ECMs(n_outlier, sorted_samples, ...
        W_including_outlier, x0_including_outlier, W_excluding_outlier, x0_excluding_outlier)
    %% This function plots the ECMs constructed including and excluding the outliers.

    % % Code written by Song Bai, Daming Li and Zhan Kang

    % For three-dimensional case, the ECMs and their projections to the coordinate planes are plotted.
    % For high-dimensional case (four and above), only the projections to the coordinate planes are plotted.

    %% Get the samples number and the dimension of the problem
    [n_samples, samples_dimension] = size(sorted_samples);

    %% Plot the ECMs
    switch samples_dimension
        case 2
            % % Plot the 2D ECMs
            figure('Name', '2D ECMs including and excluding the outliers', 'NumberTitle', 'off');

            scatter(...
                sorted_samples(end - n_outlier + 1:end, 1), ...
                sorted_samples(end - n_outlier + 1:end, 2), ...
                50, 'r', 'filled');
            hold on;

            Ellipse_plot(W_including_outlier, x0_including_outlier, 'r', 0, 200);

            scatter(...
                sorted_samples(1:n_samples - n_outlier, 1), ...
                sorted_samples(1:n_samples - n_outlier, 2), ...
                50, 'b', 'filled');

            Ellipse_plot(W_excluding_outlier, x0_excluding_outlier, 'b', 0, 200);

            plot_legend = legend('Outliers', 'ECM including outliers', 'Normal samples', 'ECM excluding outliers');
            set(plot_legend, 'Fontname', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 10);

            hold off;

        case 3
            % % Plot the 3D ECMs and their projections to the coordinate planes
            sorted_samples = sorted_samples.';

            figure('Name', '3D with projections to the coordinate planes');

            scatter3(...
                sorted_samples(1, end - n_outlier + 1:end), ...
                sorted_samples(2, end - n_outlier + 1:end), ...
                sorted_samples(3, end - n_outlier + 1:end), ...
                125, 'x', 'MarkerEdgeColor', 'r', 'LineWidth', 3);

            hold on;
            axis equal;

            %             Ellipse_plot(W_including_outlier, x0_including_outlier, 'r', 0.2);

            scatter3(...
                sorted_samples(1, 1:n_samples - n_outlier), ...
                sorted_samples(2, 1:n_samples - n_outlier), ...
                sorted_samples(3, 1:n_samples - n_outlier), ...
                55, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'b');

            Ellipse_plot(W_excluding_outlier, x0_excluding_outlier, 'b', 0.2);

            color = ["cyan", "green", "magenta"];
            flag_isnormal = true;
            plot_3D_to_2D_projection(W_excluding_outlier, x0_excluding_outlier, sorted_samples(:, 1:n_samples - n_outlier), color, flag_isnormal);

            flag_isnormal = false;
            plot_3D_to_2D_projection(W_including_outlier, x0_including_outlier, sorted_samples(:, end - n_outlier + 1:end), color, flag_isnormal);

            hold off;

            %% Plot the projections on the 2D planes separately
            x1x2_index = nchoosek(1:3, 2);
            n_figures = size(x1x2_index, 1);

            sorted_samples_transposed = sorted_samples.';

            for ii = 1:n_figures

                % % Subtract the 2D components of the samples
                sorted_samples_2D = sorted_samples_transposed(:, x1x2_index(ii, :));

                % % Get the characteristic matrix of the projection ellipse
                W_projection_including_outlier = get_w_projection(W_including_outlier, x1x2_index(ii, :));
                W_projection_excluding_outlier = get_w_projection(W_excluding_outlier, x1x2_index(ii, :));

                % % Create figure window
                figure_title = strcat("Projection on the ", num2str(x1x2_index(ii, 1)), '-', num2str(x1x2_index(ii, 2)), " plane");
                figure('Name', figure_title);

                % % Plot
                Ellipse_plot(W_projection_including_outlier, x0_including_outlier(x1x2_index(ii, :)), 'r', 0, 200);

                hold on;

                Ellipse_plot(W_projection_excluding_outlier, x0_excluding_outlier(x1x2_index(ii, :)), 'b', 0, 200);

                scatter(...
                    sorted_samples_2D(1:n_samples - n_outlier, 1), ...
                    sorted_samples_2D(1:n_samples - n_outlier, 2), ...
                    50, 'b', 'filled');

                scatter(...
                    sorted_samples_2D(end - n_outlier + 1:end, 1), ...
                    sorted_samples_2D(end - n_outlier + 1:end, 2), ...
                    50, 'r', 'filled');

                plot_legend = legend('ECM including outliers', 'Outliers', 'ECM excluding outliers', 'Normal samples');
                set(plot_legend, 'Fontname', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 10);

                hold off;
            end

        otherwise
            % % Plot the 2D projections of the high-dimensional ECMs
            x1x2_index = nchoosek(1:samples_dimension, 2);
            n_figures = size(x1x2_index, 1);

            for ii = 1:n_figures

                % % Subtract the 2D components of the samples
                sorted_samples_2D = sorted_samples(:, x1x2_index(ii, :));

                % % Get the characteristic matrix of the projection ellipse
                W_projection_including_outlier = get_w_projection(W_including_outlier, x1x2_index(ii, :));
                W_projection_excluding_outlier = get_w_projection(W_excluding_outlier, x1x2_index(ii, :));

                % % Create figure window
                figure_title = strcat("Projection on the ", num2str(x1x2_index(ii, 1)), '-', num2str(x1x2_index(ii, 2)), " plane");
                figure('Name', figure_title);

                % % Plot
                Ellipse_plot(W_projection_including_outlier, x0_including_outlier(x1x2_index(ii, :)), 'r', 0, 200);

                hold on;

                Ellipse_plot(W_projection_excluding_outlier, x0_excluding_outlier(x1x2_index(ii, :)), 'b', 0, 200);

                scatter(...
                    sorted_samples_2D(end - n_outlier + 1:end, 1), ...
                    sorted_samples_2D(end - n_outlier + 1:end, 2), ...
                    50, 'r', 'filled');

                scatter(...
                    sorted_samples_2D(1:n_samples - n_outlier, 1), ...
                    sorted_samples_2D(1:n_samples - n_outlier, 2), ...
                    50, 'b', 'filled');

                plot_legend = legend('ECM including outliers', 'Outliers', 'ECM excluding outliers', 'Normal samples');
                set(plot_legend, 'Fontname', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 10);

                hold off;
            end

    end

end
