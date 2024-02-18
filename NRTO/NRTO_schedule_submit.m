function NRBTO_schedule_submit(user_specified_time)
    %% This function submits the batch job of the NRBTO at the scheduled time

    % +-------------------------------------------------+--------------------------------------------------+
    % |                                               Input                                                |
    % +-------------------------------------------------+--------------------------------------------------+
    % | user_specified_time(optional, 22:45 by default) | Scheduled execution time (in the foramt of HHMM) |
    % +-------------------------------------------------+--------------------------------------------------+

    %%
    if nargin > 0
        user_specified_time_hour = floor(user_specified_time / 100);
        user_specified_time_minute = floor(rem(user_specified_time / 100, 1) * 100);

    else
        % % The default scheduled time is 22:45
        user_specified_time_hour = 22;
        user_specified_time_minute = 45;

    end

    user_specified_time_datenum = (user_specified_time_hour * 3600 + user_specified_time_minute * 60) / (24 * 3600);

    now_datenum = rem(now, 1);

    if user_specified_time_datenum < now_datenum
        start_delay = 0;
    else
        start_delay = floor((user_specified_time_datenum - now_datenum) * 24 * 3600);
    end

    %% Submit the batch job at the specified time
    timer_object = create_timer_callback(start_delay);
    start(timer_object);

end

% % Local functions
function timer_object = create_timer_callback(start_delay)

    timer_object = timer;
    timer_object.StartFcn = @start_function;
    timer_object.TimerFcn = @timer_function;
    timer_object.StopFcn = @stop_function;
    timer_object.ErrorFcn = @error_function;
    timer_object.StartDelay = start_delay;
end

function start_function(~, ~)

    type prompt_text.txt;

end

function timer_function(~, ~)
    NRBTO_batch_job();

end

function stop_function(mTimer, ~)
    delete(mTimer);
end

function error_function(~, ~)
    fprintf('\nSomething is wrong. See the error message please.\n');
end
