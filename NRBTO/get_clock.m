function clock_string = get_clock()
    %% This function returns the date and time as a single string
    %
    %% Input:
    %   Void
    %% Output:
    %   eg: clock_string = '202103041711'

    %%
    system_clock = clock;
    clock_string = '';

    for ii = 1:5

        [~, clock_length] = size(num2str(system_clock(ii)));

        if clock_length < 2
            clock_string_temp = strcat('0', num2str(system_clock(ii)));
        else
            clock_string_temp = num2str(system_clock(ii));
        end

        clock_string = strcat(clock_string, clock_string_temp);
    end

end
