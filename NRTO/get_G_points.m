function [G_points, G_weights] = get_G_points(G_level)
    %% This function calculates the Gaussian quadrature points and the corresponding weights

    switch G_level
        case 3
            % Three points
            Gaussian_points_1D = [-0.774596669241483, 0, 0.774596669241483];
            Gaussian_weight_1D = [0.555555555555556, 0.888888888888889, 0.555555555555556];

        case 4
            % Four points
            Gaussian_points_1D = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053];
            Gaussian_weight_1D = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454];

    end

    G_points = setprod(Gaussian_points_1D, Gaussian_points_1D);
    G_weights = kron(Gaussian_weight_1D, Gaussian_weight_1D);

end
