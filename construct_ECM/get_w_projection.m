function w_projection = get_w_projection(w, index_target)
%% This function calculates the 2D projection ellipses' characteristic matrices of the high-dimensional ellipsoid convex model

    w_swap = row_column_swap(w, index_target);

    J = w_swap(1:2, 1:2);
    L = w_swap(3:end, 1:2);
    K = w_swap(3:end, 3:end);

    w_projection = J - L.' * (K \ L);

end
