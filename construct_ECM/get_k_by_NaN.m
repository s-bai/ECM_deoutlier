function [k, NN_n, NN_unsorted, NN_D, x] = get_k_by_NaN()
    %% This function calculates the k-value using the Natural Neighborhood algorithm proposed by J Huang et al.
    % J Huang et al. A non-parameter outlier detection algorithm based on Natural Neighbor. Knowledge-based Systems. 2016.
    % https://dx.doi.org/10.1016/j.knosys.2015.10.014

    % Code written by Song Bai, Daming Li and Zhan Kang

    %     Output                                     Description
    %  ------------- ----------------------------------------------------------------------------
    %   k             The k value
    %   NN_n          The sample number of the k-nearest neighborhood
    %   NN_unsorted   The unsorted k-nearest neighbors
    %   NN_D          The relative distances between the i-th sample and its k-nearest neighbors
    %   x             The normalized samples

    % Note that the MATLAB in-house function 'knnsearch' returns the 'k numbers' of nearest neighbors,
    % which is different from the concept of the k-distance neighborhood in the present paper.
    % Thus, the in-house function 'rangesearch' is adopted for correction.

    %% Load setting parameters from the JSON file
    JSON_data = jsondecode(fileread('find_outlier_settings.json'));

    %% Load the samples data from the Excel file
    samples_table = readtable(JSON_data.samples_filename{:});
    samples = samples_table{:, :};

    %% Normalize the samples
    % The training set is stored ROW-wisely, where the i-th row contains the i-th sample
    x = mapminmax(samples.', 0, 1);
    x = x.';
    clear samples;

    % The number of the samples
    n = size(x, 1);

    % The initial k-th nearest neighbor's sequence number
    k = 1;

    % The k number for the in-house KNN search function
    k_KNN = zeros(n, 1000);
    k_KNN(:, 1) = 2;

    % The initial number of empty RN sets;
    n_empty_RN_old = 0;

    % % The times that the i-th sample is contained by the neighborhood of other samples
    n_NN = zeros(n, 1000, 'int16');

    % The sorted k-neighborhood
    NN = zeros(1000, n, 'int16');

    % The unsorted k-neighborhood
    NN_unsorted = zeros(1000, n, 'int16');

    % The relative distances of the k-nearest neighbors
    NN_D = zeros(1000, n);

    % The numbers of the k-nearest neighbors: NN_n(i, j) corresponds to the i-nearest neighbor numbers of the j-th sample
    NN_n = zeros(1000, n, 'int16');

    % The k-reverse-neighborhood
    RN = zeros(1000, n, 'int16');
    RN(1, :) = 1:n;

    % The reverse-nearest neighbors number
    RN_n = zeros(1000, n, 'int16');

    % The counter that counts the unchange number of times
    counter_unchange = 0;

    counter_unchange_threshold = 3;

    %% Build the kd tree
    Mdl = KDTreeSearcher(x);

    %% Start the loop
    while counter_unchange < counter_unchange_threshold

        for ii = 1:n

            % % Calculate the k_KNN nearest neighbors
            [~, D] = knnsearch(Mdl, x(ii, :), 'k', k_KNN(ii, k), 'IncludeTies', true);

            D = cell2mat(D);

            % % Correct the results in the presence of duplicated samples
            [neighbor_index, corrected_D] = rangesearch(x, x(ii, :), D(end) + eps);

            neighbor_index = cell2mat(neighbor_index);
            D = cell2mat(corrected_D);

            % % Calculate the corrected number of k-nearest neighbors (including duplicates)
            true_neighbor_number = length(D);
            k_KNN(ii, k) = true_neighbor_number;

            % % Store the k-nearest neighborhood
            NN_unsorted(1:true_neighbor_number, ii) = neighbor_index.';
            NN_D(1:true_neighbor_number, ii) = D.';

            NN(2:true_neighbor_number, ii) = sort(setdiff(neighbor_index, ii));
            NN(1, ii) = ii;
            NN_n(k, ii) = true_neighbor_number;

            % % Add the n_NN of the samples in NN by 1
            n_NN(setdiff(NN(1:k_KNN(ii, k), ii), ii), ii) = n_NN(setdiff(NN(1:k_KNN(ii, k), ii), ii), ii) + 1;

        end

        % % Calculate the reverse neighborhood
        for ii = 1:n

            logical_match = ismember(NN(2:end, :), ii);

            [~, match_column] = find(logical_match);
            match_column = int16(match_column);

            if ~isempty(match_column)

                reverse_neighbor_number = length(match_column);

                RN_n(k, ii) = reverse_neighbor_number;

                RN(2:reverse_neighbor_number + 1, ii) = match_column;

            end

        end

        % % Calculate the number of samples whose n_NN = 0
        n_empty_RN = n - length(find(sum(n_NN, 2)));

        % % Calculate the termination criterion
        if k == 1

            counter_unchange = counter_unchange + 1;
            n_empty_RN_old = n_empty_RN;

        else

            if n_empty_RN == n_empty_RN_old

                counter_unchange = counter_unchange + 1;

            else

                n_empty_RN_old = n_empty_RN;
                counter_unchange = 1;

            end

        end

        if counter_unchange < counter_unchange_threshold

            if k < 1000
                k_KNN(:, k + 1) = k_KNN(:, k) + 1;

            end

            %             if k<20

            k = k + 1;

            %             else
            %                return;
            %             end

        end

    end

end
