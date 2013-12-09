% Rewrite the Run file using the CellData classes
% Purpose:
%   create a simple and intuitive interface for test driving new algorithms
%   The gory details of working with any of the fcs data should be
%   abstracted away

% Formalize this as the control script
clear all; clc;

%%%%% VARIABLE INITIALIZATION %%%%%

% Cell Categories - definitions per the biologists
StemCells = Set({'HSC', 'MPP', 'CMP', 'GMP', 'MEP'});
BCells = Set({'Plasma cell', 'Pre-B I', 'Pre-B II', 'Immature B', 'Mature CD38lo B', 'Mature CD38mid B'});
TCells = Set({'Mature CD4+ T', 'Mature CD8+ T', 'Naive CD4+ T', 'Naive CD8+ T'});
NK = Set({'NK'});
pDC = Set({'Plasmacytoid DC'});
Monocytes = Set({'CD11b- Monocyte', 'CD11bhi Monocyte', 'CD11bmid Monocyte'});

% User Variables
whichCellTypes = Monocytes; 
numRandTrainExPerFile = 800; 
plotIn3D = false;
hueSensitivity = .75;
whichStimLevels = Set({'Basal'}); % Either 'Basal' or 'PV04', can contain both
useSurfaceProteinsOnly = true;

% Script Variables
cell_data_array = []; % array of pointers to all CellData objects created 

%%%%% DATA PARSING %%%%%

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
    st = cell_data_array(i).cell_stimulation_levels; % also a string
    if(~whichCellTypes.contains(ct) || ~whichStimLevels.contains(st))
        removeIndicies = [removeIndicies, i]; %#ok<AGROW>
    end
end
cell_data_array(removeIndicies) = [];

% Create single CellData object out of desired data
DesiredCells = CellData.merge(cell_data_array);

% Obtain data matrix and pre-process with arcsinh
if(useSurfaceProteinsOnly)
    data_stack = DesiredCells.getSurfaceProteinData();
else
    data_stack = DesiredCells.data;
end
data_stack = asinh(data_stack/5);

% Get data chunk indices - indices to the chunks of data that form the
% child object DesiredCells. Used for plotting
chunk_indices = DesiredCells.data_subset_indicies_for_merged_data;

% Get colors - may change for different algorithms
colors = zeros(length(whichCellTypes), 3); % RGB for every cell subtype
for j = 1:whichCellTypes.length()
    colors(j,:) = CellSubtype2Hue(whichCellTypes.list{j}, hueSensitivity);
end


%%%%% ALGORITHM SELECTION %%%%%

% Run PCA on the data
[coeff,score_PCA,latent] = princomp(data_stack);


%%%%% DATA PLOTTING %%%%%

figure;
hold on;
for i = 1:whichCellTypes.length()
    lb = chunk_indices(i);
    ub = chunk_indices(i+1)-1;
    if(plotIn3D)
        display('plotting 3D');
        scatter3(score_PCA(lb:ub,1),score_PCA(lb:ub,2),score_PCA(lb:ub,3), 20, colors(i,:));
    else
        % Plot scatter for PCA
%         subplot(1,2,1);
        scatter(score_PCA(lb:ub,1),score_PCA(lb:ub,2), 20, colors(i,:));
        title('PCA')
        
%         % Plot scatter for tSNE
%         subplot(1,2,2);
%         scatter(score_tSNE(lb:ub,1),score_tSNE(lb:ub,2), 20, colors(i,:));
%         title('t-SNE')
    end
end
legend(whichCellTypes.list)