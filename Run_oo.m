% Rewrite the Run file using the CellData classes
% Purpose:
%   create a simple and intuitive interface for test driving new algorithms
%   The gory details of working with any of the fcs data should be
%   abstracted away

% Formalize this as the control script
clear all; clc;

%% %%%%% DATA PRE-PROCESSING %%%%%

% Script Variables
cell_data_array_all = []; % array of pointers to all CellData objects created

% Create array of all fcs files as CellData objects
D=dir;
fnames = {D.name};
for i = 1:length(fnames)
    if(CellData.isReadable(fnames{i}))
        cell_data_array_all = [cell_data_array_all, CellData(fnames{i})]; %#ok<AGROW>
    end
end

%% DATA PROCESSING

%%%%% VARIABLE INITIALIZATION %%%%%

% Cell Categories - definitions per the biologists
StemCells = Set({'HSC', 'MPP', 'CMP', 'GMP', 'MEP'});
BCells = Set({'Plasma cell', 'Pre-B I', 'Pre-B II', 'Immature B', 'Mature CD38lo B', 'Mature CD38mid B'});
TCells = Set({'Mature CD4+ T', 'Mature CD8+ T', 'Naive CD4+ T', 'Naive CD8+ T'});
NK = Set({'NK'});
pDC = Set({'Plasmacytoid DC'});
Monocytes = Set({'CD11b- Monocyte', 'CD11bhi Monocyte', 'CD11bmid Monocyte'});

% User Variables
whichCellTypes = Monocytes & pDC & NK & TCells & BCells; 
% whichCellTypes = TCells;
numRandTrainExPerFile = 200; %seems optimal for tsne 
hueSensitivity = 2.5;
whichStimLevels = Set({'Basal'}); % Either 'Basal' or 'PV04', can contain both
useSurfaceProteinsOnly = true;


%%%%% DATA PARSING %%%%%

% Keep the CellData objects whose cell_type is contained in the
% whichCellTypes variable
removeIndicies = [];
for i = 1:length(cell_data_array_all)
    ct = cell_data_array_all(i).cell_types; % will be a single string since no CellData objects have seen merger
    st = cell_data_array_all(i).cell_stimulation_levels; % also a string
    if(~whichCellTypes.contains(ct) || ~whichStimLevels.contains(st))
        removeIndicies = [removeIndicies, i]; %#ok<AGROW>
    end
end
cell_data_array = cell_data_array_all;
cell_data_array(removeIndicies) = [];

% Create single CellData object out of desired data
DesiredCells = CellData.merge(cell_data_array, numRandTrainExPerFile);

% Obtain data matrix and pre-process with arcsinh
if(useSurfaceProteinsOnly)
    data_stack = DesiredCells.getSurfaceProteinData();
else
    data_stack = DesiredCells.getProteinData();
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

%% t-SNE %%
% dimensionality reduction to dim = 2
% see the file alg_tsne for more details

% Want dimensionality reduction to 2
dim = 2;

% stopping criteria: number of iterations is no more than 100, runtime is
% no more than 30 seconds, and the relative tolerance in the embedding is 
% no less than 1e-3. Taken from Max's tsne example demo_swissroll.m
opts.maxit = 400; opts.runtime = 900; opts.tol = 1e-3;
opts.X0 = 1e-5*randn(size(data_stack, 1), dim);

% Run algorithm
[tsne_output, E, A, T] = alg_tsne(data_stack, dim, opts);

% Plot results
figure; hold on;
for i = 1:whichCellTypes.length()
    lb = chunk_indices(i);
    ub = chunk_indices(i+1)-1;        
    scatter(tsne_output(lb:ub,1),tsne_output(lb:ub,2), 20, colors(i,:));
    title(['TSNE: iter #' num2str(length(E)), ', e=' num2str(E(end)),...
       ', t=' num2str(T(end)), ', N/file=' num2str(numRandTrainExPerFile)]);   
end
legend(whichCellTypes.list)

%% s-SNE %%
% dimensionality reduction to dim = 2
% see the file alg_ssne for more details

% Want dimensionality reduction to 2
dim = 2;

% stopping criteria: number of iterations is no more than 100, runtime is
% no more than 30 seconds, and the relative tolerance in the embedding is 
% no less than 1e-3. Taken from Max's tsne example demo_swissroll.m
opts.maxit = 600; opts.runtime = 900; opts.tol = 1e-3;
opts.X0 = 1e-5*randn(size(data_stack, 1), dim);

% Run algorithm
[ssne_output, E, A, T] = alg_ssne(data_stack, dim, opts);

% Plot results
figure; hold on;
for i = 1:whichCellTypes.length()
    lb = chunk_indices(i);
    ub = chunk_indices(i+1)-1;        
    scatter(ssne_output(lb:ub,1),ssne_output(lb:ub,2), 20, colors(i,:));
    title(['SSNE: iter #' num2str(length(E)), ', e=' num2str(E(end)),...
       ', t=' num2str(T(end)), ', N/file=' num2str(numRandTrainExPerFile)]);   
end
legend(whichCellTypes.list)

%% EE %%
% dimensionality reduction to dim = 2
% see the file alg_ee for more details

% Want dimensionality reduction to 2
dim = 2;

% stopping criteria: number of iterations is no more than 100, runtime is
% no more than 30 seconds, and the relative tolerance in the embedding is 
% no less than 1e-3. Taken from Max's tsne example demo_swissroll.m
opts.maxit = 100; opts.runtime = 900; opts.tol = 1e-3;
opts.X0 = 1e-5*randn(size(data_stack, 1), dim);

% Run algorithm
[ee_output, E, A, T] = alg_ee(data_stack, dim, opts);

% Plot results
figure; hold on;
for i = 1:whichCellTypes.length()
    lb = chunk_indices(i);
    ub = chunk_indices(i+1)-1;        
    scatter(ee_output(lb:ub,1),ee_output(lb:ub,2), 20, colors(i,:));
    title(['EE: iter #' num2str(length(E)), ', e=' num2str(E(end)),...
       ', t=' num2str(T(end)), ', N/file=' num2str(numRandTrainExPerFile)]);   
end
legend(whichCellTypes.list)

%% PCA %%
[coeff,score_PCA,latent] = princomp(data_stack);
plotIn3D = false;
figure; hold on;
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
        title(['PCA: N/file=' num2str(numRandTrainExPerFile)]);
        
%         % Plot scatter for tSNE
%         subplot(1,2,2);
%         scatter(score_tSNE(lb:ub,1),score_tSNE(lb:ub,2), 20, colors(i,:));
%         title('t-SNE')
    end
end
legend(whichCellTypes.list)
