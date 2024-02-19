function make_movie(pause_interval_user_specified, record_flag)
    %% This function makes movie of topological evolution from the latest result file

    % % Input
    % pause_interval_user_specified: Pause interval between two frames (0.05 sec by default)
    % record_flag: The flag marking whether to make the movie (1 for true, 0 for false)

    %% Load the latest NRTO results MAT file
    mat_file_infor = dir('NRTO*.mat');

    mat_file_timestamp = zeros(length(mat_file_infor), 1);

    for ii = 1:length(mat_file_infor)

        mat_file_timestamp(ii) = datenum(mat_file_infor(ii).date);

    end

    [~, mat_file_timestamp_latest] = max(mat_file_timestamp);

    fprintf('\n\n\nThe latest NRTO results file is: %s\n', mat_file_infor(mat_file_timestamp_latest).name);
    results = load(mat_file_infor(mat_file_timestamp_latest).name);

    %%
    max_iter = results.iter;

    if nargin > 0
        pause_interval = pause_interval_user_specified;
    else
        pause_interval = 0.05;
    end

    figure;
    set(gca, 'nextplot', 'replacechildren');

    if nargin > 1

        if record_flag ~= 0
            video_handle = VideoWriter('NRTO_iteration_movie.avi');
            video_handle.FrameRate = 24;
            open(video_handle);

        end

    end

    for ii = 1:max_iter
        colormap(gray);
        imagesc(flipud(1 - results.rho_history(:, :, ii)));
        caxis([0 1]);
        axis equal;
        axis off;
        drawnow;

        iter_string = num2str(ii);
        annotation_text = strcat('Iteration #', iter_string);

        text(0, -15, annotation_text);

        pause(pause_interval);

        if nargin > 1

            if record_flag ~= 0
                frame = getframe(gcf);
                writeVideo(video_handle, frame);
            end

        end

    end

    if nargin > 1

        if record_flag ~= 0
            close(video_handle);
        end

    end

end
