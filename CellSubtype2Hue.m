function rgb = CellSubtype2Hue(cellSubtype, sensitivity)

typeToColorMap = {'StemCell', 'red';...
                    'B', 'green';...
                    'T', 'blue';...
                    'NK', 'black';...
                    'pDC', 'yellow';...
                    'Monocyte', 'magenta'
    };

cellType = CellTypeFromSubCell(cellSubtype);

for i = 1:size(typeToColorMap, 1)
    if(strcmpi(cellType, typeToColorMap{i,1}))
        rgb = GenerateRandomColorHue(typeToColorMap{i,2}, sensitivity);
    end
end
