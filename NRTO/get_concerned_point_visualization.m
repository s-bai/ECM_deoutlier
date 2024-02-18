function new_plot_handle_concerned_point = get_concerned_point_visualization(figure_number_CPts, iter_index, plot_handle_concerned_point, concerned_point, n_ECMs, Ui)
    %% This function returns the visualization of the concerned points on the standard ellipsoid convex models
    %
    %% ---------------------------------------------- Input ----------------------------------------------
    % figure_number_CPts: The figure numbers
    % iter_index: The iteration step
    % plot_handle_concerned_point: The graphic object handles
    % concerned_point: The concerned points of the current iteration
    % n_ECMs: The number of the ellipsoid convex model
    % Ui:   The displacement constraints and the corresponding reliability constraints
    % n_Ui: The number of the concerned displacements
    %
    %% ---------------------------------------------- Output ----------------------------------------------
    % new_plot_handle_concerned_point: The updated graphic object handles

    %%
    [n_Ui, ~] = size(Ui);

    eta_lower_limit = Ui(:, 3);

    %% Plot the background standard ellipsoid convex models in the first iteration step
    if iter_index == 1

        for ii = 1:n_Ui
            figure(figure_number_CPts(ii));

            x = 0; y = 0; r = eta_lower_limit(ii);

            theta = 0:pi / 50:2 * pi;
            xunit = r * cos(theta) + x;
            yunit = r * sin(theta) + y;
            
            plot(xunit, yunit, 'LineWidth', 2, 'Color', [0.3010, 0.7450, 0.9330]);
            
            xlim([-0.5 - eta_lower_limit(ii), 1 + eta_lower_limit(ii)]);
            ylim([-0.5 - eta_lower_limit(ii), 1 + eta_lower_limit(ii)]);
            
            pbaspect([1 1 1]);

        end

    end

    %%
    new_plot_handle_concerned_point = cell(n_ECMs, n_Ui);

    marker_color = [1, 0, 0; 0, 1, 0; 0, 0, 1];

    legend_content = ["ECM*", "CPt#1", "CPt#2", "CPt#3", "CPt#4"];

    for ii = 1:n_Ui

        figure(figure_number_CPts(ii));
        hold on

        flag_delete_former_scatter = 1;

        for jj = 1:n_ECMs

            % Erase the former concerned points on the graphics after the 1st iteration
            if iter_index > 1

                if flag_delete_former_scatter

                    for kk = 1:n_ECMs
                        delete(plot_handle_concerned_point{kk, ii});
                    end

                    flag_delete_former_scatter = 0;
                end

            end

            % % Plot the current concerned points on the standard ECM
            new_plot_handle_concerned_point{jj, ii} = scatter(...
                concerned_point(2 * jj - 1, ii), concerned_point(2 * jj, ii), 40, marker_color(jj, :), 'filled');
            legend(legend_content(1:jj + 1));
        end

    end

end
