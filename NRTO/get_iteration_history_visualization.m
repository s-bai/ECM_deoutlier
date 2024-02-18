function get_iteration_history_visualization(figure_number_iter_history, iter, iteration_history_array, y_limits)
    %% This function plots the iteration history dynamically
    %
    %% ---------------------------------------------- Input ----------------------------------------------
    % figure_number_iter_history: The handle number of the graphic object
    % iter: The iteration step
    % iteration_history_array: The vector storing the iteration history
    % y_limits: The lower and upper limits of the y-axis
    %
    %% ---------------------------------------------- Output ----------------------------------------------
    % Void

    %% Switch to the target figure window
    figure(figure_number_iter_history);
    hold on;

    %% Plot 
    plot(1:iter, iteration_history_array(1:iter), 'LineWidth', 2, 'Color', 'b');
    xlim([0, iter + 10]);
    ylim([0.7 * y_limits(1), 1.1 * y_limits(2)]);

end
