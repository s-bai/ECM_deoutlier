function n_outlier = detect_outlier_MAD(LOF_values)
    %% This function detect the outliers by identifying the samples with abnormally large LOF scores.
    % The Median Absolute Deviation measurement is adopted.

    % % Code written by Song Bai, Daming Li and Zhan Kang

    %% Variables initialization
    ii = length(LOF_values);
    n_outlier = 0;
    exit_flag = false;

    %% Identify all samples with abnormally LOF scores. 
    %  Only the abnormally large LOF scores indicate outliers.
    %  The in-house function isoutlier of MATLAB is employed.
    outlier_logical_tag = isoutlier(LOF_values);

    %% 
    while ~exit_flag

        if outlier_logical_tag(ii)
            n_outlier = n_outlier + 1;

            ii = ii - 1;

        else
            exit_flag = true;
        end

        if ii == 0
            exit_flag = true;
        end

    end

end
