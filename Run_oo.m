% Rewrite the Run file using the CellData classes
% Purpose:
%   create a simple and intuitive interface for test driving new algorithms
%   The gory details of working with any of the fcs data should be
%   abstracted away

% Formalize this as the control script
clear all; clc;

% Cell Categories - definitions per the biologists
StemCells = Set({'HSC', 'MPP', 'CMP', 'GMP', 'MEP'});
BCells = Set({'Plasma cell', 'Pre-B I', 'Pre-B II', 'Immature B', 'Mature CD38lo B', 'Mature CD38mid B'});
TCells = Set({'Mature CD4+ T', 'Mature CD8+ T', 'Naive CD4+ T', 'Naive CD8+ T'});
NK = Set({'NK'});
pDC = Set({'Plasmacytoid DC'});
Monocytes = Set({'CD11b- Monocyte', 'CD11bhi Monocyte', 'CD11bmid Monocyte'});

% User Variables
whichCellTypes = pDC; 
numRandTrainExPerFile = 800; 
plotIn3D = false;
hueSensitivity = 5;
cell_stimulation_level = 'Basal'; % Either 'Basal' or 'PV04'

% Script Variables
cell_data_array = []; % array of pointers to all CellData objects created 

% Create array of all fcs files as CellData objects
D=dir;
fnames = {D.name};
for i = 1:length(fnames)
    if(CellData.isReadable(fnames{i}))
        cell_data_array = [cell_data_array, CellData(fnames{i})]; %#ok<AGROW>
    end
end

% Keep the CellData objects whose cell_type is contained in the
% whichCellTypes variable
removeIndicies = [];
for i = 1:length(cell_data_array)
    ct = cell_data_array(i).cell_types; % will be a single string since no CellData objects have seen merger
    if(~whichCellTypes.contains(ct))
        removeIndicies = [removeIndicies, i]; %#ok<AGROW>
    end
end
cell_data_array(removeIndicies) = []