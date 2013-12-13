function mergeFiguresIntoSubplot(nrows, ncols, figHandles)

% Get all open figure handles
figHandles_all = findobj('Type','figure'); % Get all

% Get subplot positions
nfindex = max(figHandles_all) + 1;
nf = figure(nfindex); % Create new figure
sppos = []
for i = 1:length(figHandles)
    ax = subplot(nrows, ncols,i);
    sppos = [sppos; get(ax, 'pos')];
    cla(ax);
end
cla(nf);

% Copy figures into subplots of new figure
new_splots = {};
for i = 1:length(figHandles)
    new_splots{end +1} = copyobj(get(figHandles(i), 'children'), nfindex);
end
for i = 1:length(figHandles)
    set(new_splots{i}, 'pos', sppos(i,:));
end