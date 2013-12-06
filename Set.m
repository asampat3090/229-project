
% Cell array that acts like a mathematical set
% no duplicates allowed
classdef Set < handle
    
    properties
        list = {};
    end
    
    methods
        function obj = Set(items)
            if(~iscell(items))
                items = {items};
            end
            obj.list = items;
        end
        
        % Returns true if this contains the item
        % if the item is a cell, it returns true if this contains every
        % element of that cell
        function tf = contains(this, item)
            tf = false;
            if(iscell(item))
                tf = true;
                for i = 1:length(item)
                    tf = tf & this.contains(item{i});
                end
            else
                for i = 1:length(this.list)
                    if(isequal(this.list{i}, item))
                        tf = true;
                    end
                end
            end
        end
        
        function add(this, item)
            if(iscell(item))
                for i = 1:length(item)
                    this.add(item{i})
                end
            else
                if(~this.contains(item))
                    this.list{end+1} = item;
                end                
            end
            
        end
%         
%         % make this recursive WHENCE
%         function add(this, items)
%             if(~iscell(items))
%                 items = {items};
%             end
%             
%             for i = 1:length(items)
%                 if(~this.contains(items{i}))
%                     this.list{end+1} = items{i};
%                 end
%                 
%             end
%         end
    end
    
end