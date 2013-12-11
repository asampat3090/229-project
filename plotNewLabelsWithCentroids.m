function plotNewLabelsWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, plot_title, new_labels)

chunk_indices = DesiredCells.data_subset_indicies_for_merged_data;

% PLOT - New Labels with Centroids
% display('Plotting PCA - New Labels')
figure('name',figName); hold on;
for i = 1:length(DesiredCells.cell_types)
    lb = chunk_indices(i);
    ub = chunk_indices(i+1)-1;
    for j = lb:ub
        % Plot scatter for PCA
        ct = CellData.possible_cell_types{new_labels(j)};        
        r_idx = DesiredCells.RelativeIndexOfCellType(ct);
        rgb = DesiredCells.colors_relative(r_idx,:);
        scatter(processed_data(j,1),processed_data(j,2), 20, rgb, 'fill');
    end
end
title(plot_title);
legend(DesiredCells.cell_types)
for i = 1:size(centroids, 1) % Centroids
    if(centroids(i,1) ~= inf)
        ct = CellData.possible_cell_types{i};
        r_idx = DesiredCells.RelativeIndexOfCellType(ct);
        rgb = DesiredCells.colors_relative(r_idx,:);
        scatter(centroids(i,1), centroids(i,2), 400, rgb, 'h', 'fill', 'markeredgecolor', abs(1-rgb));
    end
end
drawnow