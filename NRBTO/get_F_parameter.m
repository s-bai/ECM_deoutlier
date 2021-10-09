function F_parameter = get_F_parameter(nelx, nely, inv_T, F_DOFs_and_F0)
    %% This function calculates the sparse inverse of all transformation matrices and the sparse nominal loads vector
    %
    %% Input:
    %       nelx, nely: Numbers of elements in x and y directions
    %
    %       inv_T{3}: The inverse of the transformation matrices of the uncertain loads
    %
    %       F_DOFs_and_F0: The DOFs (the 1st column) and the nominal values (the 2nd column) of applied loads
    %
    %% Output:
    %       F_parameter.inv_T_sparse 
    %
    %       F_parameter.F_nominal_sparse 
    
    %%
    [~, n_inv_T] = size(inv_T);
    n_inv_T_F = n_inv_T - 2;

    %% Assign the inverse matrices of the transformation matrices of the uncertain loads to inv_T_F
    inv_T_F = cell(1, n_inv_T_F);
    F_nominal = cell(1, n_inv_T_F);

    for ii = 1:n_inv_T_F
        inv_T_F{ii} = inv_T{ii + 2};
        F_nominal{ii} = F_DOFs_and_F0(2 * ii - 1:2 * ii, 2);
    end

    F_parameter.inv_T = zeros(4 * n_inv_T_F, 1);
    F_parameter.F_DOF = zeros(2 * n_inv_T_F, 1);
    F_parameter.F_nominal = zeros(2 * n_inv_T_F, 1);

    F_parameter.inv_T(1:4) = reshape(inv_T_F{1}, [], 1);
    F_parameter.F_nominal(1:2, 1) = F_nominal{1};

    for ii = 1:n_inv_T_F
        F_parameter.F_DOF(2 * ii - 1:2 * ii, 1) = F_DOFs_and_F0(2 * ii - 1:2 * ii, 1);
    end

    for ii = 2:n_inv_T_F
        F_parameter.inv_T(4 * ii - 3:4 * ii) = reshape(inv_T_F{ii}, [], 1);
        F_parameter.F_nominal(1:2 * ii, 1) = cat(1, F_parameter.F_nominal(1:2 * (ii - 1), 1), F_nominal{ii});
    end

    % Create the sparse matrix which consists of the elements of T inverse
    temp_DOF = reshape(F_parameter.F_DOF.', 2, n_inv_T_F);
    temp_DOF = kron(temp_DOF, ones(1, n_inv_T_F));

    iii = reshape(temp_DOF, 1, []);
    jjj = kron(1:2 * n_inv_T_F, ones(1, n_inv_T_F));

    F_parameter.inv_T_sparse = sparse(iii, jjj, F_parameter.inv_T, 2 * (nelx + 1) * (nely + 1), 2 * n_inv_T_F);
    F_parameter.F_nominal_sparse = sparse(F_parameter.F_DOF, ones(size(F_parameter.F_DOF)), F_parameter.F_nominal, 2 * (nelx + 1) * (nely + 1), 1);

end
