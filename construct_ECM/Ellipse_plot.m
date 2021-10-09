function Ellipse_plot(A, C, color, ECM_opacity, varargin)
    %% This function plots 2D ellipse or 3D ellipsoid

    % % This function is modified based on the open-source code 'Ellipse_plot' developed by Nima Moshtagh

    % % The license for the current code is attached at the bottom of this file

    % % Modified by Song Bai

    %  Ellipse_Plot(A,C,N) plots a 2D ellipse or a 3D ellipsoid
    %  represented in the "center" form:
    %
    %                   (x-C)' A (x-C) <= 1
    %
    %  A and C could be the outputs of the function: "MinVolEllipse.m",
    %  which computes the minimum volume enclosing ellipsoid containing a
    %  set of points in space.
    %
    %  Inputs:
    %  A: a 2x2 or 3x3 matrix.
    %  C: a 2D or a 3D vector which represents the center of the ellipsoid.
    %  color: The line color of the 2D ellipse
    %  ECM_opacity: The opacity of the 3D ellipsoid
    %  N: the number of grid points for plotting the ellipse; Default: N = 20.
    %
    %  Nima Moshtagh
    %  nima@seas.upenn.edu
    %  University of Pennsylvania
    %  Feb 1, 2007
    %  Updated: Feb 3, 2007

    %%%%%%%%%%%  start  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = 20; % Default value for grid

    % See if the user wants a different value for N.
    %----------------------------------------------
    if nargin > 4
        N = varargin{1};
    end

    % check the dimension of the inputs: 2D or 3D
    %--------------------------------------------
    if length(C) == 3
        Type = '3D';
    elseif length(C) == 2
        Type = '2D';
    else
        disp('Cannot plot an ellipse with more than 3 dimensions!!');
        return
    end

    % "singular value decomposition" to extract the orientation and the
    % axes of the ellipsoid
    [~, D, V] = svd(A);

    if strcmp(Type, '2D')
        % get the major and minor axes
        %------------------------------------
        a = 1 / sqrt(D(1, 1));
        b = 1 / sqrt(D(2, 2));

        theta = 0:1 / N:2 * pi + 1 / N;

        % Parametric equation of the ellipse
        %----------------------------------------
        state(1, :) = a * cos(theta);
        state(2, :) = b * sin(theta);

        % Coordinate transform
        %----------------------------------------
        X = V * state;
        X(1, :) = X(1, :) + C(1);
        X(2, :) = X(2, :) + C(2);

    elseif strcmp(Type, '3D')
        % generate the ellipsoid at (0,0,0)
        %----------------------------------
        a = 1 / sqrt(D(1, 1));
        b = 1 / sqrt(D(2, 2));
        c = 1 / sqrt(D(3, 3));
        [X, Y, Z] = ellipsoid(0, 0, 0, a, b, c, N);

        %  rotate and center the ellipsoid to the actual center point
        %------------------------------------------------------------
        XX = zeros(N + 1, N + 1);
        YY = zeros(N + 1, N + 1);
        ZZ = zeros(N + 1, N + 1);

        for k = 1:length(X)

            for j = 1:length(X)
                point = [X(k, j) Y(k, j) Z(k, j)]';
                P = V * point;
                XX(k, j) = P(1) + C(1);
                YY(k, j) = P(2) + C(2);
                ZZ(k, j) = P(3) + C(3);
            end

        end

    end

    % Plot the ellipse
    %----------------------------------------
    if strcmp(Type, '2D')
        line_specification = strcat('-', color);
        plot(X(1, :), X(2, :), line_specification, 'LineWidth', 1.5);
        hold on;

        axis equal
    else
        surface_color = zeros(N + 1, N + 1, 3);

        switch color
            case 'r'
                surface_color(:, :, 1) = ones(N + 1);
            case 'b'
                surface_color(:, :, 3) = ones(N + 1);
        end

        ECM_opacity = num2str(ECM_opacity);

        ellipsoid_surface = surf(XX, YY, ZZ, surface_color, 'FaceAlpha', ECM_opacity);
        ellipsoid_surface.EdgeColor = 'flat';
        % ellipsoid_surface.FaceColor = 'flat';
        ellipsoid_surface.LineStyle = 'none';
        ellipsoid_surface.FaceLighting = 'gouraud';

        box on;
        set(gca, 'Projection', 'perspective', 'BoxStyle', 'full');

        axis equal
        hidden off
    end

    %%
    %     Copyright (c) 2009, Nima Moshtagh
    % All rights reserved.

    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are
    % met:

    %     * Redistributions of source code must retain the above copyright
    %       notice, this list of conditions and the following disclaimer.
    %     * Redistributions in binary form must reproduce the above copyright
    %       notice, this list of conditions and the following disclaimer in
    %       the documentation and/or other materials provided with the distribution

    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.
