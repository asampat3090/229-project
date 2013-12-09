% File type: script
% Purpose: 
    % Parse all fcs files in the directory and gather list of all the cell
    % types in those files. 
    % requires that in the file name the part following
    % the last '_' but before the '.' of the extension is the cell type name
    % concatenates all these cell types into a cell variable and ensures no
    % repeats.
% Returns:
%   cell variable 'CellTypes' - a list of strings, each string is a
%   biological cell type contained in the directory.

function CellTypes = cellTypesFromDirectory();
    D=dir;
    fnames = {D.name};

    CellTypes = {};
    for i = 1:length(fnames)

        % Take file in directory, obtain the filename without extension, and if
        % it isn't a .fcs file, skip this file
        file = fnames{i};
        [file, ext] = strtok(file, '.');
        if(~strcmpi(ext, '.fcs'))
            continue;
        end

        % Scan through the filename to get the last part (after the last '_' ),
        % which contains the cell type name
        [str, rem] = strtok(file, '_');
        while(~isempty(rem))
            [str, rem] = strtok(rem, '_');        
        end
        celltype = str;

        % Check that CellTypes doesn't already contain the new cell type
        isNewCell = true;
        for j = 1:length(CellTypes)
            if(strcmpi(celltype, CellTypes{j}))
                isNewCell = false;
            end
        end

        % IF its a new cell type, add it to the list
        if(isNewCell)
%             display(['Adding Cell ' celltype]);
            CellTypes{end+1} = celltype;
        end

    end

end