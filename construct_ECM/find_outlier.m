function find_outlier()
    %% This is the main function of the present research project
    %  This function detects the outliers and plots the ECMs

    % % Code written by Song Bai, Daming Li and Zhan Kang

    %% Prepare for results output
    results.values = cell(8, 1);

    %% Load settings from the JSON file
    JSON_data = jsondecode(fileread('find_outlier_settings.json'));

    %% Load the samples data from the Excel file
    samples_table = readtable(JSON_data.samples_filename{:});
    samples = samples_table{:, :};

    % The flag for normalization
    normalization = JSON_data.normalization;

    %% Calculate the k value by the NaN method
    [k, NN_n, NN_unsorted, NN_D, normalized_samples] = get_k_by_NaN();

    fprintf('\n\n--------------------------------------------------------------\n\n');
    fprintf('The k is obtained by the NaN algorithm. The k value is k = %d', k);
    fprintf('\n\n--------------------------------------------------------------\n\n');

    %% Get the LOF values
    LOF_values_and_samples = get_LOF_value(samples, normalized_samples, k, NN_n, NN_unsorted, NN_D, normalization);

    %% Detect the outliers by the MAD method
    n_outlier = detect_outlier_MAD(LOF_values_and_samples(:, 2));

    %% Calculate the W and x0 of the ECMs
    [n_outlier, LOF_values_and_samples, w_including_outlier, x0_including_outlier, w_excluding_outlier, x0_excluding_outlier] = ...
        get_w_x0(LOF_values_and_samples, n_outlier);

    %% Plot the ECMs
    plot_ECMs(n_outlier, LOF_values_and_samples(:, 3:end), w_including_outlier, x0_including_outlier, w_excluding_outlier, x0_excluding_outlier);

    % % Get the upper and lower limits of the samples including/excluding the outliers
    limits_including_outliers = [...
                                min(samples); ...
                                max(samples) ...
                                ];

    limits_excluding_outliers = [...
                                min(samples(LOF_values_and_samples(1:length(LOF_values_and_samples(:, 1)) - n_outlier, 1), :)); ...
                                max(samples(LOF_values_and_samples(1:length(LOF_values_and_samples(:, 1)) - n_outlier, 1), :)) ...
                                ];

    % % Calculate the volume ratio of the ECMs constructed using the samples including & excluding the outliers
    A = cell(1, 2);
    w = cell(1, 2);
    w{1} = w_including_outlier;
    w{2} = w_excluding_outlier;

    for ii = 1:2
        [v_temp, d_temp] = eig(w{ii});
        A{ii} = v_temp * sqrt(d_temp) * v_temp.';
    end

    volume_ratio_gamma = det(A{2}) / det(A{1});

    %% Organize the results
    results.values{1} = w_including_outlier;
    results.values{2} = x0_including_outlier;
    results.values{3} = limits_including_outliers(1, :);
    results.values{4} = limits_including_outliers(2, :);

    results.values{5} = w_excluding_outlier;
    results.values{6} = x0_excluding_outlier;
    results.values{7} = limits_excluding_outliers(1, :);
    results.values{8} = limits_excluding_outliers(2, :);

    results.values{9} = volume_ratio_gamma;

    results.titles = [...
                    "The matrix W of the ECM including the outliers is:", ...
                    "The x0 of the ECM including the outliers is:", ...
                    "The lower limit of the samples including the ourliers is:", ...
                    "The upper limit of the samples including the ourliers is:", ...
                    "The matrix W of the ECM excluding the outliers is:", ...
                    "The x0 of the ECM excluding the outliers is:", ...
                    "The lower limit of the samples excluding the ourliers is:", ...
                    "The upper limit of the samples excluding the ourliers is:", ...
                    "The volume ratio is:"
                ];

    %% Print the results to screen
    for ii = 1:9

        fprintf('\n%s\n', results.titles(ii));
        disp(results.values{ii});

    end

    fprintf('\nThe number of outliers is n = %d\n', n_outlier);

    %% Save the results to the file
    filename = JSON_data.samples_filename{:};
    filename = strcat('result_', filename(1:end - 4), 'mat');

    save(filename, ...
        'w_including_outlier', 'x0_including_outlier', 'limits_including_outliers', ...
        'w_excluding_outlier', 'x0_excluding_outlier', 'limits_excluding_outliers', ...
        'volume_ratio_gamma', 'n_outlier', 'k');

end
