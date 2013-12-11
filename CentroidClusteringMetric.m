function [label_err new_labels centroids] = CentroidClusteringMetric(metric_data, DesiredCells)

% CLUSTERING METRIC BASED ON AVERAGING
% metric_data = score_PCA(:,1:2);

%Initialize centroids in 2D
% initialize centroid for all cell subtypes
% if that subtype isn't considered in this run, its centroid is set to inf
%  -- Makes it easier to bookkeep which centroid corresponds to which cell
%  type
% the ith row of mu will correspond to the centroid of the cell subtype:
%   CellData.possible_cell_types{i};

chunk_indices = DesiredCells.data_subset_indicies_for_merged_data;
mu = Inf*ones(length(CellData.possible_cell_types), size(metric_data,2)); 
for i = 1:length(DesiredCells.cell_types)
    lb = chunk_indices(i); ub = chunk_indices(i+1)-1;
    data_chunk = metric_data(lb:ub,:);
    idx = CellData.AbsoluteIndexOfCellType(DesiredCells.cell_types{i});
    mu(idx,:) = sum(data_chunk,1) ./ size(data_chunk,1);
end

% Generate New Labels that reference the CellData.possible_cell_types list
data_distances = pdist2(metric_data, mu, 'euclidean'); % matrix of distance: (i,j) is distance from metric_data(i,:) to mu(j,:)
[nonsense new_labels] = min(data_distances, [], 2); %new_labels is column vector
new_labels = new_labels';

% Compare to Actual Labels
actual_labels = DesiredCells.data_celltype_indices;
label_err = sum(actual_labels ~= new_labels)/length(actual_labels);

centroids = mu;