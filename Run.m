clear all;
clc;
currentFolder = pwd;
addpath(genpath(currentFolder));

% Tutorial:
% set the MarrowRegion, set whichCellTypes, set numRandTrainExPerFile, set
% proteinRange, 

% for i =1:length(CellTypes)
%     display(['hue of ' CellTypes{i} ' is ' num2str(CellSubtype2Hue(CellTypes{i}))])
% end

%% Get Cell Types from all Basal Files in directory

% User variables
% MarrowRegion = 'PVO4';
MarrowRegion = 'Basal';

D=dir(['*' MarrowRegion '*']);
fnames = {D.name};
% CellTypes = {'Pre-B II', 'Pre-B I', 'Plasmacytoid DC', 'Plasma cell', 'NK', 'Naive CD8+ T',...
%                 'blue ball'
%     };

% Grab all Cell Types - requires that in the file name the part followed
% the last '_' but before the '.' of the extension is the cell type name
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
        display(['Adding Cell ' celltype]);
        CellTypes{end+1} = celltype;
    end
    
end
%% Read and Plot data for desired cell types

% Define the cell types of interest & initialize vector of the number of
% data pts per cell used (will be useful in color coding the PCA plots)

% Cell Super-Types (note, I switch between supercelltype and cell type all the
% time...fix this in a later code version
StemCells = {'HSC', 'MPP', 'CMP', 'GMP', 'MEP'};
BCells = {'Plasma cell', 'Pre-B I', 'Pre-B II', 'Immature B', 'Mature CD38lo B', 'Mature CD38mid B'};
TCells = {'Mature CD4+ T', 'Mature CD8+ T', 'Naive CD4+ T', 'Naive CD8+ T'};
NK = {'NK'};
pDC = {'Plasmacytoid DC'};
Monocytes = {'CD11b- Monocyte', 'CD11bhi Monocyte', 'CD11bmid Monocyte'};

% User Variables
% whichCellTypes = CellTypes; 
whichCellTypes = Monocytes; 
numRandTrainExPerFile = 800; 
plotIn3D = false;
hueSensitivity = 5;

proteinRange = 3:40; %Should be within 3:40 (3 is the first column index corresponding to a protein, 40 is the last). With Karen I did 5:30
useSurfaceProteins = false; %if set to true, the script will only consider proteins with 'CD' in their names

% Gather data from files which correspond to the cell types of interest
numCellsPerCellType = zeros(1,length(whichCellTypes));
dataStack = [];
colors = zeros(length(whichCellTypes), 3);
for j = 1:length(whichCellTypes)
    wc = whichCellTypes{j};
    
    %Generate a color hue for this subcell type based on its general cell
    %type
    colors(j,:) = CellSubtype2Hue(wc, hueSensitivity);
    
    % If this file is of a cell type of interest, add it to the datastack
    % and update the number of data pts per cell being used
    
    for i=1:length(fnames)
        f = fnames{i};
        if(strfind(f, wc))
%             display(['Cell Type ' wc ': Using file ' f]);
            [data, header] = fca_readfcs(f); 
            
            % Determine the range as numRandTrainExPerFile random integers
            % between 1 and the number of cell data points in this file
            numtrainExamples = size(data,1);
            display(num2str(numtrainExamples));
            r = randperm(numtrainExamples);
            range = r(1:min(numRandTrainExPerFile, numtrainExamples));
            
            % Use Surface Proteins if flagged
            if(useSurfaceProteins)
                proteinNames = {header.par.name2};
                surfaceProteinIndices = [];
                for k = 1:length(proteinNames)
                    if(~isempty(findstr(proteinNames{k}, 'CD')))
                        surfaceProteinIndices = [surfaceProteinIndices k];
                    end
                end
                proteinRange = surfaceProteinIndices; 
            end
            
            dataStack = [dataStack; data(range,proteinRange)];
            numCellsPerCellType(j) = numCellsPerCellType(j) + length(range);
%             range = randint(1,400, [1 size(data,1)]);
%             dataStack = [dataStack; data(range,5:30)];
        end            
    end
end

% Run PCA on the data
dataStack = asinh(dataStack/5);
[coeff,score_PCA,latent] = princomp(dataStack);
size(coeff)
size(score_PCA)
size(dataStack)
size(latent)

% Run tSNE on the data  - still unsure of output

% The function performs symmetric t-SNE on the NxD dataset X to reduce its 
% dimensionality to no_dims dimensions (default = 2). The data is 
% preprocessed using PCA, reducing the dimensionality to initial_dims 
% dimensions (default = 30). Alternatively, an initial solution obtained 
% from an other dimensionality reduction technique may be specified in 
% initial_solution. The perplexity of the Gaussian kernel that is employed 
% can be specified through perplexity (default = 30). The labels of the
% data are not used by t-SNE itself, however, they are used to color
% intermediate plots. Please provide an empty labels matrix [] if you
% don't want to plot results during the optimization.
% The low-dimensional data representation is returned in mappedX.

score_tSNE = tsne(dataStack);
%TSNE Performs symmetric t-SNE on dataset X
%
%   mappedX = tsne(X, labels, no_dims, initial_dims, perplexity)
%   mappedX = tsne(X, labels, initial_solution, perplexity)
%

% Initialize score indices based on numCellsPerCellType (for group
% scattering)
scoreIndices = zeros(1, length(numCellsPerCellType));
for i = 1:length(scoreIndices)
    scoreIndices(i) = sum(numCellsPerCellType(1:i));
end
scoreIndices = [0, scoreIndices];

% Create figure and scatter plot the data, labelling colors based on cell
% types

% Add each scatter plot as a subplot in grid.

figure;
hold on;
for i = 1:length(whichCellTypes)
    lb = scoreIndices(i)+1;
    ub = scoreIndices(i+1);
    if(plotIn3D)
        display('plotting 3D');
        scatter3(score_PCA(lb:ub,1),score_PCA(lb:ub,2),score_PCA(lb:ub,3), 20, colors(i,:));
    else
        % Plot scatter for PCA
        subplot(1,2,1);
        scatter(score_PCA(lb:ub,1),score_PCA(lb:ub,2), 20, colors(i,:));
        title('PCA')
        % Plot scatter for tSNE
        subplot(1,2,2);
        scatter(score_tSNE(lb:ub,1),score_tSNE(lb:ub,2), 20, colors(i,:));
        title('t-SNE')
    end
%     scatter(score((i-1)*400+1:400*i,1),score(i:400*i + 1,2),'filled','b');
%     scatter(score((i-1)*400+1:400*i,1),score((i-1)*400+1:400*i,2), c(i),'filled');
%     scatter(score(lb:ub,1),score(lb:ub,2), c(i),'filled', 'markersize', 10);
end
legend(whichCellTypes)



space = ' ';
newline = '\n'; % for use with sprintf
if(useSurfaceProteins)
    str_sp = ' surface proteins only';
else 
    str_sp = ' all proteins';
end
expr_title = sprintf([MarrowRegion space 'data' str_sp newline '<' num2str(numRandTrainExPerFile) ' random cells taken from each file']);
xlabel('p1');
ylabel('p2');
title(expr_title);
set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
set(gca, 'fontsize', 14, 'linewidth', 2);

