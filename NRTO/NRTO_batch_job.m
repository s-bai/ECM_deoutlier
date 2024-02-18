function NRTO_batch_job()
    %% This script submits the batch job of the NRTO numerical examples (void input argument)
    NRTO_batch_filenames = dir('NRTO_cantilever_beam_*.json');
    n_files = length(NRTO_batch_filenames);

    %% Get the current directory
    FolderName = pwd;

    %% Submit the batch job
    for ii = 1:n_files

        nrto_lof_2d(NRTO_batch_filenames(ii).name);

        %% Save all opened figures as .fig file
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

        for iFig = 1:length(FigList)
            FigHandle = FigList(iFig);
            FigName = get(FigHandle, 'Name');

            savefig(FigHandle, fullfile(FolderName, strcat(FigName, '_job#', num2str(ii), '.fig')));
        end

        %% Close all opened figures
        close all;
    end

end
