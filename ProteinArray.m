classdef ProteinArray < handle
    %ProteinARray A smart array object that contains all the proteins used and
    %their corresponding tag metals. Has the ability to reference either
    %field by the other
    %   Detailed explanation goes here
    
    properties
%         Protein.name;
%         Protein.metal;
        len = 0;
    end
    
    methods
        function addPair(this, p, m)
            if(~ischar(p) || ~ischar(m))
                error('Protein Array pairs must be strings');
            end
            if(this.
            this.proteins{this.len+1} = p;
            this.metals{this.len+1} = m;
            this.len = this.len + 1;
        end
        
        % Returns true of the protein array contains the protein p
        function tf = containsProtein(this, p)
            tf = false;
            for i = 1:this.len
                if(strcmpi(p, this.proteins{i}))
                    tf = true;
                    break;
                end
            end
        end
        
        % Returns true of the protein array contains the metal m
        function tf = containsMetal(this, m)
            tf = false;
            for i = 1:this.len
                if(strcmpi(m, this.metals{i}))
                    tf = true;
                    break;
                end
            end
        end 
    end
    
end



% classdef ProteinArray < handle
%     %ProteinARray A smart array object that contains all the proteins used and
%     %their corresponding tag metals. Has the ability to reference either
%     %field by the other
%     %   Detailed explanation goes here
%     
%     properties
%         proteins = {};
%         metals = {};
%         len = 0;
%     end
%     
%     methods
%         function addPair(this, p, m)
%             if(~ischar(p) || ~ischar(m))
%                 error('Protein Array pairs must be strings');
%             end
%             if(this.
%             this.proteins{this.len+1} = p;
%             this.metals{this.len+1} = m;
%             this.len = this.len + 1;
%         end
%         
%         % Returns true of the protein array contains the protein p
%         function tf = containsProtein(this, p)
%             tf = false;
%             for i = 1:this.len
%                 if(strcmpi(p, this.proteins{i}))
%                     tf = true;
%                     break;
%                 end
%             end
%         end
%         
%         % Returns true of the protein array contains the metal m
%         function tf = containsMetal(this, m)
%             tf = false;
%             for i = 1:this.len
%                 if(strcmpi(m, this.metals{i}))
%                     tf = true;
%                     break;
%                 end
%             end
%         end 
%     end
%     
% end
% 
