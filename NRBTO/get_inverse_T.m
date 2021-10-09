function [n_T, inv_T] = get_inverse_T(w)
    %% This function calculates the inverse matrices of all characteristic matrices
    %
    %% Input:
    %   w{1}: The charateristic matrix of the uncertain length and width
    %   w{2}: The charateristic matrix of the uncertain Young's modulus and the Poisson's ratio
    %   w{3}: The charateristic matrix of the uncertain loads
    %
    %% Output:
    %   n_T: Number of the ellipsoid convex models
    %   inv_T{1}: The inverse of the charateristic matrix of the uncertain length and width
    %   inv_T{2}: The inverse of the charateristic matrix of the uncertain Young's modulus and the Poisson's ratio
    %   inv_T{3}: The inverses of the charateristic matrix of the uncertain loads

    %%
    [~, n_T] = size(w);
    inv_T = cell(size(w));

    for ii = 1:n_T
        [eigenvector_w, eigenvalue_w] = eig(w{ii});

        inv_T{ii} = inv(eigenvector_w * sqrt(eigenvalue_w) * eigenvector_w.');
    end

end
