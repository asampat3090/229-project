%<<<<<<< HEAD
c = CellData('Marrow1_15_PVO4_Marrow1_PVO4_Mature CD8+ T.fcs');
d = CellData('Marrow1_15_PVO4_Marrow1_PVO4_GMP.fcs');
e = CellData('Marrow1_01_Basal1_Marrow1_Basal1_Naive CD8+ T.fcs');
m_cd = CellData.merge([c d]);
m_de = CellData.merge([d e]);
m = CellData.merge([m_cd m_de]);

%%
clear all; clc;
cell = {1 2 3}
s = Set(cell)
s.contains(1)
s.contains(2)
s.contains({1})
s.contains({1, {2}})
currentFolder = pwd;
addpath(genpath(currentFolder));

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
hueSensitivity = 10;

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

%% 
%%%%%%%%%%%%%% Prepare the Data for Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Setup all of the relevant strings 
space = ' ';
newline = '\n'; % for use with sprintf
if(useSurfaceProteins)
    str_sp = ' surface proteins only';
else 
    str_sp = ' all proteins';
end
expr_title = sprintf([MarrowRegion space 'data' str_sp newline '<' num2str(numRandTrainExPerFile) ' random cells taken from each file']);


%%%%%%%%% DIMENSIONAL REDUCTION ALGORITHMS & PLOTTING %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Naive Linear Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% Run PCA on the data
display('Running PCA on dataset')
dataStack = asinh(dataStack/5);
[coeff,score_PCA,latent] = princomp(dataStack);
size(coeff)
size(score_PCA)
size(dataStack)
size(latent)

% Plot the Figure for PCA
display('Plotting PCA Result')
figure('name','PCA');
hold on;
for i = 1:length(whichCellTypes)
    lb = scoreIndices(i)+1;
    ub = scoreIndices(i+1);
    if(plotIn3D)
        display('plotting 3D');
        scatter3(score_PCA(lb:ub,1),score_PCA(lb:ub,2),score_PCA(lb:ub,3), 20, colors(i,:));
    else
        % Plot scatter for PCA
        scatter(score_PCA(lb:ub,1),score_PCA(lb:ub,2), 20, colors(i,:));
        xlabel('p1');
        ylabel('p2');
        title(strcat('PCA ',expr_title));
    end
%     scatter(score((i-1)*400+1:400*i,1),score(i:400*i + 1,2),'filled','b');
%     scatter(score((i-1)*400+1:400*i,1),score((i-1)*400+1:400*i,2), c(i),'filled');
%     scatter(score(lb:ub,1),score(lb:ub,2), c(i),'filled', 'markersize', 10);
end
legend(whichCellTypes)
hold off;
drawnow
%% 

% % Run non-classical MDS on the data
% display ('Running non-classical MDS')
% dissimilarities = pdist(zscore(dataStack));
% [score_ncMDS,stress] = mdscale(dissimilarities,2,'criterion','metricstress');
% 
% % Plot the Figure for non-classical MDS
% display ('Plotting non-classical MDS Result')
% figure('name','ncMDS');
% hold on;
% for i = 1:length(whichCellTypes)
%     lb = scoreIndices(i)+1;
%     ub = scoreIndices(i+1);
%     if(plotIn3D)
%         display('plotting 3D');
%         scatter3(score_ncMDS(lb:ub,1),score_ncMDS(lb:ub,2),score_ncMDS(lb:ub,3), 20, colors(i,:));
%     else
%         % Plot scatter for PCA
%         scatter(score_ncMDS(lb:ub,1),score_ncMDS(lb:ub,2), 20, colors(i,:));
%         xlabel('p1');
%         ylabel('p2');
%         title(strcat('non-classical MDS ',expr_title));
%     end
% %     scatter(score((i-1)*400+1:400*i,1),score(i:400*i + 1,2),'filled','b');
% %     scatter(score((i-1)*400+1:400*i,1),score((i-1)*400+1:400*i,2), c(i),'filled');
% %     scatter(score(lb:ub,1),score(lb:ub,2), c(i),'filled', 'markersize', 10);
% end
% legend(whichCellTypes)
% hold off;
% drawnow
%% 

% Run ICA on the data - good for separation

%%%%%%%%%%%%%% Non - Linear Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% Run Isomap Algorithm on Data
display('Running Isomap on Data');
[score_isomap, mapping_isomap] = isomap(dataStack);
%ISOMAP Runs the Isomap algorithm
%
%   [mappedX, mapping] = isomap(X, no_dims, k); 

% Plot the Figure for Isomap
display('Plotting Isomap Result')
figure('name','Isomap');
hold on;
for i = 1:length(whichCellTypes)
    lb = scoreIndices(i)+1;
    ub = scoreIndices(i+1);
    if ub>size(score_isomap,1)
        ub=size(score_isomap,1)
    end
    if(plotIn3D)
        display('plotting 3D');
        if(size(score_isomap,2)<3)
            display('Cannot plot SNE results in 3D - need more data - Run tSNE with no_dims of >=3');
        else
            scatter3(score_isomap(lb:ub,1),score_isomap(lb:ub,2),score_isomap(lb:ub,3), 20, colors(i,:));
        end
    else
        % Plot scatter for tSNE
        scatter(score_isomap(lb:ub,1),score_isomap(lb:ub,2), 20, colors(i,:));
        xlabel('p1');
        ylabel('p2');
        title(strcat('Isomap ',expr_title));
    end
%     scatter(score((i-1)*400+1:400*i,1),score(i:400*i + 1,2),'filled','b');
%     scatter(score((i-1)*400+1:400*i,1),score((i-1)*400+1:400*i,2), c(i),'filled');
%     scatter(score(lb:ub,1),score(lb:ub,2), c(i),'filled', 'markersize', 10);
end
legend(whichCellTypes)
hold off;
drawnow

%%%%%%%%%%% SNE & t-SNE ALGORITHMS (take a while to converge) %%%%%%%%%%%
%% 

% % Run SNE on Data 
% display('Running SNE on Data')
% score_SNE = sne(dataStack);
% %SNE Implementation of Stochastic Neighbor Embedding
% %
% %   mappedX = sne(X, no_dims, perplexity)

% % Plot the Figure for SNE
% display('Plotting SNE Result')
% figure('name','SNE');
% hold on;
% for i = 1:length(whichCellTypes)
%     lb = scoreIndices(i)+1;
%     ub = scoreIndices(i+1);
%     if(plotIn3D)
%         display('plotting 3D');
%         if(size(score_SNE,2)<3)
%             display('Cannot plot SNE results in 3D - need more data - Run tSNE with no_dims of >=3');
%         else
%             scatter3(score_SNE(lb:ub,1),score_SNE(lb:ub,2),score_SNE(lb:ub,3), 20, colors(i,:));
%         end
%     else
%         % Plot scatter for tSNE
%         scatter(score_SNE(lb:ub,1),score_SNE(lb:ub,2), 20, colors(i,:));
%         xlabel('p1');
%         ylabel('p2');
%         title(strcat('SNE ',expr_title));
%     end
% %     scatter(score((i-1)*400+1:400*i,1),score(i:400*i + 1,2),'filled','b');
% %     scatter(score((i-1)*400+1:400*i,1),score((i-1)*400+1:400*i,2), c(i),'filled');
% %     scatter(score(lb:ub,1),score(lb:ub,2), c(i),'filled', 'markersize', 10);
% end
% legend(whichCellTypes)
% hold off;
% drawnow
%% 

% % Run tSNE on the data 
% display('Running t-SNE on data')
% score_tSNE = tsne(dataStack);
% %TSNE Performs symmetric t-SNE on dataset X
% %
% %   mappedX = tsne(X, labels, no_dims, initial_dims, perplexity)
% %   mappedX = tsne(X, labels, initial_solution, perplexity)

% % Plot the Figure for t-SNE
% display('Plotting t-SNE Result')
% figure('name','t-SNE');
% hold on;
% for i = 1:length(whichCellTypes)
%     lb = scoreIndices(i)+1;
%     ub = scoreIndices(i+1);
%     if(plotIn3D)
%         display('plotting 3D');
%         if(size(score_tSNE,2)<3)
%             display('Cannot plot tSNE results in 3D - need more data - Run tSNE with no_dims of >=3');
%         else
%             scatter3(score_tSNE(lb:ub,1),score_tSNE(lb:ub,2),score_tSNE(lb:ub,3), 20, colors(i,:));
%         end
%     else
%         % Plot scatter for tSNE
%         scatter(score_tSNE(lb:ub,1),score_tSNE(lb:ub,2), 20, colors(i,:));
%         xlabel('p1');
%         ylabel('p2');
%         title(strcat('t-SNE ',expr_title));
%     end
% %     scatter(score((i-1)*400+1:400*i,1),score(i:400*i + 1,2),'filled','b');
% %     scatter(score((i-1)*400+1:400*i,1),score((i-1)*400+1:400*i,2), c(i),'filled');
% %     scatter(score(lb:ub,1),score(lb:ub,2), c(i),'filled', 'markersize', 10);
% end
% legend(whichCellTypes)
% hold off;
% drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Run k-means on UNSUPERVISED data AFTER reduced to 2 or 3 dimensions
% score_kmeans = kmeans(dataStack,3);
% size(score_kmeans)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
set(gca, 'fontsize', 14, 'linewidth', 2);

