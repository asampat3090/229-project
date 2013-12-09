% Celldata
% this class reads a .fcs file and stores its data in a matrix
% Purpose: abstract away the logistical handling of the fcs data and be
%   able to merge multiple sets of data as needed
%
% Notation Convention: classes and object instances are capitalized
%                       all other variables are strictly lower case and use
%                       '_' delineators
%                       functions are initially lower case, but if they are
%                       the union of multiple words, then the 2nd+ words
%                       are capitalized
%
% FCS file data format
%   contains a large matrix where rows correspond to individual cells and
%   the columns correspond to some variable. This variable is usually a
%   protein, in which case the matrix entry would be the counts of that
%   protein in the processed cell. Other variables include time, etc.
%
%   contains a header struct, which labels the columns of the matrix as
%   well as contains data like the name of the cell type

classdef CellData < handle
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    properties(Constant = true)
        fcs_reader = @fca_readfcs; % object designed to work for this particular reader
        fcs_ext = '.fcs';
        possible_stimulation_levels = {'Basal', 'PVO4'}; 
        surface_protein_tag = 'CD'; % case sensitive
    end
    
    properties
        data = []; % [m x n] matrix read from the .fcs file        
        column_headings = {}; % string of the names of the columns
        cell_stimulation_levels = {}; % Basal, PV04
        cell_types = {}; % HSC, GMP, MPP, Mature CD38mid B, etc.
    end
    
    properties(SetAccess = private)
        data_subset_indicies_for_merged_data = [1]; % if obj formed under merging, it contains row indices to the points at which the data changes from one CellData parent to the nex
    end
    
    properties(Access = private)
        protein_column_index = 3; % the fcs data has columns 3:end as protein data, the first two columns are Time and Cell Length
        header = {}; % data struct as read in by the reader
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% STATIC CLASS METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    methods(Static = true)
        
        % ISREADABLE %
        % Checks that 'filename' is a string with:
        % suffix = '.fcs'
        % 
        function tf = isReadable(filename)
            tf = true;
            if(~ischar(filename))
                tf = false;
                return;
            end
            
            [~, ext] = strtok(filename, '.');
            if(~strcmpi(ext, CellData.fcs_ext))
                tf = false;
                return
            end            
        end
        
        % MERGE %
        % Creates a child object from all the parents in pointer_array
        % pointer_array is an array of pointers to ClassData objects
        %
        % throws an error if the column_headings are not identical (case
        % sensitive)
        % the child has:
        %   the concatenated data of the parents (in no
        %   particular order) along the 1st dimension of the matrix (cell
        %   rows)
        %   the cell_stimulation_levelss of all the parents 
        %   the cell_typess of all the parents
        % the child does NOT have:
        %   the same header as the parents (in fact, it'll have an empty
        %   header)
        % 
        % the integer array data_subset_indicies_for_merged_data contains
        % indices to the row indices that denote the change from one
        % parents data to the next. The last entry = size(child.data,1)+1;
        function obj = merge(pointer_array)       
            
            % Check the objects in the pointer_array are column_heading
            % compatible
            if(~isequal(pointer_array.column_headings))
                error('pointers being merged must have identical column headers');
            end
            
            % Create new object and add variable values
            obj = CellData();
            obj.column_headings = pointer_array(1).column_headings;  
            for i = 1:length(pointer_array)
                % Update data, stimulation level, and cell type
                obj.addData( pointer_array(i).data );                
                obj.addStimLevel( pointer_array(i).cell_stimulation_levels );
                obj.addCellType( pointer_array(i).cell_types );
                
                % Update data indicies
                m = size(pointer_array(i).data, 1);
                temp = obj.data_subset_indicies_for_merged_data;
                obj.data_subset_indicies_for_merged_data = [temp, temp(end) + m];
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% OBJECT METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %Constructor method - reads the fcs file 'fcs_filename' and parses
        %its data into the properties
        function this = CellData(fcs_filename)
            
            %Empty constructor (MATLAB does not support overloading fncs            
            if(nargin == 0)
                return;
            end
            
            % Parse data and header
            [this.data, this.header] = CellData.fcs_reader(fcs_filename);
            this.column_headings = {this.header.par.name2};
            
            % Parse cell stimulation level
            for i = 1:length(CellData.possible_stimulation_levels)
                if(strfind(lower(fcs_filename) , lower(CellData.possible_stimulation_levels{i})))
                    this.cell_stimulation_levels = CellData.possible_stimulation_levels(i);
                    break;
                end
            end
            
            % Parse cell type
            [fcs_filename, ~] = strtok(fcs_filename, '.');
            [str, rem] = strtok(fcs_filename, '_');
            while(~isempty(rem))
                [str, rem] = strtok(rem, '_');         %#ok<STTOK>
            end
            this.cell_types = {str};
        end
              
       
        
        % ADD DATA %
        % Adds the [m x n] matrix new_data to the data field of 'this'. 
        % m can take any value
        % the n of new_data must match the n of this.data
        function addData(this, new_data)
            this.data = [this.data; new_data];
        end
        
        % ADD CELL STIMULATION LEVEL %
        % adds a new cell stimulation level if its not already of that type
        % new_stims can be a single item or a cell of items        
        function addStimLevel(this, new_stim)
            set = Set(this.cell_stimulation_levels);
            set.add(new_stim);
            this.cell_stimulation_levels = set.list;
        end
        
        % ADD CELL TYPE %
        % adds a new cell stimulation level if its not already of that type
        function addCellType(this, new_type)
            set = Set(this.cell_types);
            set.add(new_type);
            this.cell_types = set.list;
        end
        
        % GET SURFACE PROTEINS DATA %
        % returns a matrix of values corresponding to the surface protein
        % data only
        function matrix = getSurfaceProteinData(this)
            proteinNames = this.column_headings;
            surfaceProteinIndices = [];
            for k = 1:length(proteinNames)
                if(~isempty(strfind(proteinNames{k}, CellData.surface_protein_tag)))
                    surfaceProteinIndices = [surfaceProteinIndices k]; %#ok<AGROW>
                end
            end
            matrix = this.data(:, surfaceProteinIndices);            
        end
    end
    
end