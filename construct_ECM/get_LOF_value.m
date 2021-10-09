function LOF_values_and_samples = get_LOF_value(samples, normalized_samples, k, NN_n, NN_unsorted, NN_D, normalization)
    %% This function calculates the LOF of all samples

    % % Code written by Song Bai, Daming Li and Zhan Kang

    %         Input                                          Description
    %  -------------------- -----------------------------------------------------------------------------
    %   samples              The original samples
    %   normalized_samples   The normalized samples
    %   k                    The k-value obtained by the NaN method
    %   NN_n                 The samples in the k-nearest neighborhoods
    %   NN_unsorted          The unsorted k-nearest neighbors
    %   NN_D                 The relative distances between the i-th samples and its k-nearest neighbors
    %   normalization        The flag for samples normalization

    %% The number of the samples
    n_samples = size(samples, 1);

    %% Variable initialization
    % The LOF scores
    LOF_values = zeros(n_samples, 1);

    % The regularization parameter epsilon
    epsilon = 1e-3;

    % The k-nearest reachability distances
    dR_k = zeros(1000, n_samples);

    % The average k-nearest reachability distances
    dR_k_average = zeros(n_samples, 1);

    % % Get the k-nearest distances of the samples
    d_k = zeros(n_samples, 1);

    for ii = 1:n_samples
        d_k(ii) = NN_D(NN_n(k, ii), ii);
    end

    % The indices of the k-nearest neighbors
    NN_index = zeros(1000, n_samples);

    %% Calcualte the k-nearest rechability distances and the average k-nearest reachability distances
    for ii = 1:n_samples

        % % Get the indices of the k-nearest neighbors of the ii-th sample
        NN_index_temp = NN_unsorted(1:NN_n(k, ii), ii);
        NN_index(1:NN_n(k, ii) - 1, ii) = setdiff(NN_index_temp, ii);

        % % Set the reference samples vector to the duplicates of the ii-th sample
        ref_sample = repmat(samples(ii, :), [NN_n(k, ii) - 1, 1]);

        % Calculate the real distances between the ii-th sample and the samples in k-NN
        real_distance = vecnorm(ref_sample - samples(NN_index(1:NN_n(k, ii) - 1, ii), :), 2, 2);

        % Calculate the k-nearest reachability distance of the ii-th sample
        dR_k(1:NN_n(k, ii) - 1, ii) = max([real_distance, d_k(NN_index(1:NN_n(k, ii) - 1, ii))], [], 2);

        % % Calculate the average k-nearest reachability distances of the ii-th samples
        dR_k_average(ii, 1) = mean(dR_k(1:NN_n(k, ii) - 1, ii));

    end

    %% Calcualte the LOF scores
    for ii = 1:n_samples

        % % Calculate the harmonic mean of the average k-nearest reachability distances of the samples in NN_k(ii)
        dR_k_harmonic_mean = harmmean(dR_k_average(NN_index(1:NN_n(k, ii) - 1, ii), 1));

        % % Calculate the LOF score of the ii-th sample
        LOF_values(ii, 1) = (epsilon + dR_k_average(ii)) / (epsilon + dR_k_harmonic_mean);

    end

    %% Normalize the samples if required
    if normalization

        output_samples = normalized_samples;

    else

        output_samples = samples;

    end

    % Concatanate the LOF values and the samples
    LOF_values_and_samples = [(1:n_samples).', LOF_values, output_samples];

    % Sort the samples ascendingly by the LOF values
    LOF_values_and_samples = sortrows(LOF_values_and_samples, 2);
end
