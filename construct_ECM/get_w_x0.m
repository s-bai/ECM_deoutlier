function [n_outlier, sorted_samples, W_including_outlier, x0_including_outlier, W_excluding_outlier, x0_excluding_outlier] = ...
        get_w_x0(sorted_samples, n_outlier)
    %% This function calculates the W and x0 of the ECMs

    % % Code written by Song Bai, Daming Li and Zhan Kang

    % % Note:
    % %   Some variables initializations for the CVX subroutine trigger false Code Analyzer alarm.
    % %   Thus '%#ok' are added to the end of the corresponding lines to surpress these false alarms.

    %       Input                               Description
    %  ---------------- ------------------------------------------------------------
    %   sorted_samples   The samples sorted ascendingly according to the LOF scores
    %   n_outlier        The number of the outliers

    %          Output                                  Description
    %  ---------------------- -------------------------------------------------------------
    %   n_outlier              The number of the outliers after correction
    %   sorted_samples         The sorted samples after removing the mis-detected outliers
    %   W_including_outlier    The W of the ellipsoid convex model including outliers
    %   x0_including_outlier   The center of the ellipsoid convex model including outliers
    %   W_excluding_outlier    The W of the ellipsoid convex model excluding outliers
    %   x0_excluding_outlier   The center of the ellipsoid convex model excluding outliers

    %% Get the samples number and the problem dimension
    [n_samples, problem_dimension] = size(sorted_samples(:, 3:end)); %#ok

    %% ------------------------ Modelling of the ECM including the outliers ------------------------
    % % Solve the SDP problem using the CVX
    cvx_begin

    variable A_including_outlier(problem_dimension, problem_dimension) symmetric
    variable b_including_outlier(problem_dimension)

    minimize(-det_rootn(A_including_outlier))

    subject to
    norms(A_including_outlier * sorted_samples(:, 3:end).' - b_including_outlier * ones(1, n_samples), 2) <= 1; %#ok

    cvx_end

    %%
    % The x0 of the ECM including the outliers
    x0_including_outlier = A_including_outlier \ b_including_outlier;

    % The w of the ECM including the outliers
    W_including_outlier = A_including_outlier.' * A_including_outlier;

    %% ------------------------ Modelling of the ECM excluding the outliers ------------------------
    false_outlier_counter = 0;

    % % Solve the SDP problem using the CVX
    cvx_begin

    variable A_excluding_outlier(problem_dimension, problem_dimension) symmetric
    variable b_excluding_outlier(problem_dimension)

    minimize(-det_rootn(A_excluding_outlier))

    subject to
    norms(A_excluding_outlier * sorted_samples(1:n_samples - n_outlier, 3:end).' - b_excluding_outlier * ones(1, n_samples - n_outlier), 2) <= 1; %#ok

    cvx_end

    %% Check if any outlier is identified by mistake
    query_outlier = sorted_samples(n_samples - n_outlier + 1:end, :);

    % The x0 of the ECM excluding the outliers
    x0_excluding_outlier = A_excluding_outlier \ b_excluding_outlier;

    % The w of the ECM excluding the outliers
    W_excluding_outlier = A_excluding_outlier.' * A_excluding_outlier;

    %%
    for ii = 1:n_outlier

        if (query_outlier(ii, 3:end) - x0_excluding_outlier.') * W_excluding_outlier * (query_outlier(ii, 3:end).' - x0_excluding_outlier) <= 1

            false_outlier_counter = false_outlier_counter + 1;

            temp_sample = query_outlier(ii, :);
            query_outlier(ii, :) = query_outlier(false_outlier_counter, :);
            query_outlier(false_outlier_counter, :) = temp_sample;

        end

    end

    sorted_samples(n_samples - n_outlier + 1:end, :) = query_outlier;
    n_outlier = n_outlier - false_outlier_counter;

end
