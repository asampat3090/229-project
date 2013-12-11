Broken Code

function plotDataWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, labels)

chunk_indices = DesiredCells.data_subset_indicies_for_merged_data;
colors = DesiredCells.getRelativeColorsByAbsoluteLabels(labels);

% PLOT
% display('Plotting PCA - True Labels')
figure('name',figName); hold on;
for i = 1:length(DesiredCells.cell_types)
    lb = chunk_indices(i);
    ub = chunk_indices(i+1)-1;

    % Plot scatter for PCA
    scatter(processed_data(lb:ub,1), processed_data(lb:ub,2), 20, colors(i,:));        
    
end
legend(DesiredCells.cell_types)
title(['PCA: N/file=' num2str(numRandTrainExPerFile)]);

if(~isempty(centroids))
    for i = 1:size(centroids, 1) % Centroids
        if(centroids(i,1) ~= inf)
            ct = CellData.possible_cell_types{i};
            r_idx = DesiredCells.RelativeIndexOfCellType(ct);
            rgb = DesiredCells.colors_relative(r_idx,:);
            scatter(centroids(i,1), centroids(i,2), 400, rgb, 'h', 'fill', 'markeredgecolor', abs(1-rgb));
        end
    end
end
drawnow