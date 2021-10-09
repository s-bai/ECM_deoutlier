function [Alpha_constraint, void_eq_constraint, D_Alpha_constraint, void_Deq_constraint] = get_Alpha_constraint(phi, index_concerned, Ui, n_ellipsoid)
    %% Calculate constraint function of the inner loop problem
    %
    %% ---------------------------------------------- Input ----------------------------------------------
    %
    %   phi: The standard uncertain parameters vector
    %       phi(1), phi(2): For the length and width uncertainties
    %       phi(3), phi(4): For the Young's modulus and the Poisson's ratio uncertainties
    %       phi(5), phi(6): For the loads uncertainties  
    %
    %   index_concerned: The index of the concerned structural performance
    %
    %            The arrangement of the elements of Ui
    %     ┌─────────────────┬────────────────────────┬─────────────────────────────┐
    %     │    DOFs of Ui   │ The upper limits of Ui │   The lower limits of eta   │
    %     ├─────────────────┼────────────────────────┼─────────────────────────────┤
    %     │     U1_DOF      │       U1_upper         │        eta_U1_lower         │
    %     │     U2_DOF      │       U2_upper         │        eta_U2_lower         │
    %     └─────────────────┴────────────────────────┴─────────────────────────────┘
    %
    %   n_ellipsoid: Number of ellipsoid convex models
    %
    %% ---------------------------------------------- Output ----------------------------------------------
    %
    %   Alpha_constraint: Value of nonlinear constraint function
    %   D_Alpha_constraint: Sensitivity of nonlinear constraint function with respect to phi
    %
    %   void_eq_constraint: Void equality constraint
    %   void_Deq_constraint: Void sensitivity of equality constraint

    %%
    [n_phi, ~] = size(phi);

    %% Calculate the alpha values
    Alpha_constraint = zeros(n_ellipsoid, 1);

    for ii = 1:n_ellipsoid
        Alpha_constraint(ii) = norm(phi(2 * ii - 1:2 * ii))^2 - Ui(index_concerned, 3)^2;
    end

    %%
    D_Alpha_constraint = zeros(n_phi, n_ellipsoid);

    for ii = 1:n_ellipsoid
        D_Alpha_constraint(2 * ii - 1:2 * ii, ii) = 2 * phi(2 * ii - 1:2 * ii);
    end

    %%
    void_eq_constraint = [];
    void_Deq_constraint = [];
end
