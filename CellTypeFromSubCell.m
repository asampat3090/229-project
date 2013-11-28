function cellType = CellTypeFromSubCell(subCellType)

types = {'StemCell', 'B', 'T', 'NK', 'pDC', 'Monocyte'};

StemCells = {'HSC', 'MPP', 'CMP', 'GMP', 'MEP'};
BCells = {'Plasma cell', 'Pre-B I', 'Pre-B II', 'Immature B', 'Mature CD38lo B', 'Mature CD38mid B'};
TCells = {'Mature CD4+ T', 'Mature CD8+ T', 'Naive CD4+ T', 'Naive CD8+ T'};
NK = {'NK'};
pDC = {'Plasmacytoid DC'};
Monocytes = {'CD11b- Monocyte', 'CD11bhi Monocyte', 'CD11bmid Monocyte'};
    
for i =1:length(StemCells)
    if(strcmpi(StemCells{i}, subCellType))
        cellType = types(1);
    end
end

for i =1:length(BCells)
    if(strcmpi(BCells{i}, subCellType))
        cellType = types(2);
    end
end

for i =1:length(TCells)
    if(strcmpi(TCells{i}, subCellType))
        cellType = types(3);
    end
end

for i =1:length(NK)
    if(strcmpi(NK{i}, subCellType))
        cellType = types(4);
    end
end

for i =1:length(pDC)
    if(strcmpi(pDC{i}, subCellType))
        cellType = types(5);
    end
end

for i =1:length(Monocytes)
    if(strcmpi(Monocytes{i}, subCellType))
        cellType = types(6);
    end
end

end