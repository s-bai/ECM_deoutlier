function FE_assembly = get_FE_assembly(nelx, nely)
    %% This function calculates the finite element assembly information for FEA

    node_numbers_matrix = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
    node_DOF_vector = reshape(2 * node_numbers_matrix(1:end - 1, 1:end - 1) + 1, nelx * nely, 1);
    node_DOF_matrix = repmat(node_DOF_vector, 1, 8) + repmat([0 1 2 * nely + [2 3 0 1] -2 -1], nelx * nely, 1);
    FE_assembly.iK = reshape(kron(node_DOF_matrix, ones(8, 1))', 64 * nelx * nely, 1);
    FE_assembly.jK = reshape(kron(node_DOF_matrix, ones(1, 8))', 64 * nelx * nely, 1);
    FE_assembly.node_DOF_matrix = node_DOF_matrix;

end
