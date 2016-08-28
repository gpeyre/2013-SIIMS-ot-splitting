%
% Constructor for staggered object using the dimension vector (or number) supplied
%
%
% CALL : staggered(dimvect [...],boundaryc 'Dirichlet',testt 0 1 or 2)
%
% EX   : staggered([3 6],'Dirichlet',2)
%
% NB   : defined fields : + obj.M{k}     : values describing the coordinate k of the vector field
%                         + obj.dim      : dimensions of the central grid
%                         + obj.dims{k}  : dimension  of the grid of coordinate k of the vector field obj.M{k}
%                         + obj.boundary : champ permettant la gestion des contraintes de bord du
%                                          champ de vecteur (.type, .ind{k} and  .valb{k} such that
%                                          obj.M{k}(obj.boundary.ind{k}) = obj.boundary.valb{k})
%

classdef staggered
    properties ( GetAccess = 'public', SetAccess = 'public' )
        dim                          = [];
        dims                         = {};
        M                            = {};
        boundary                     = struct;
    end
    
    methods
        function obj = staggered(dimvect,boundaryc,testt)
            
            if(exist('dimvect','var')==0)
                dimvect                  = [3 4];
            end
            if(exist('boundaryc','var')==0)
                boundaryc                = '';
            end
            if(exist('testt','var')==0)
                testt                    = 0;
            end
            obj.dim                    = dimvect;
            for k = 1:length(dimvect)
                lt                       = dimvect;
                lt(k)                    = dimvect(k) + 1;
                if(testt>0)
		  obj.M{k}               = 2*rand(lt)-1;
                else
		  obj.M{k}               = zeros(lt);
                end
                obj.dims{k}              = lt;
            end
            obj.boundary.type            = boundaryc;
	    if(isempty(obj.boundary.type)==0)
              % boundary indicis
	      for k = 1:length(obj.M)
                Ivt                      = ind2subv(size(obj.M{k}),1:numel(obj.M{k}));
                st                       = size(obj.M{k});
                It1                      = find(Ivt(:,k)==1);
                It2                      = find(Ivt(:,k)==st(k));
                obj.boundary.ind{k}      = [It1(:);It2(:)];
                obj.boundary.valb{k}     = zeros(numel(obj.boundary.ind{k}),1);
                if(testt>1)
                    obj.boundary.valb{k} = k;
                else
                    obj.boundary.valb{k} = 0;
                end
	      end
            end
            obj                        = projectonboundaryconditions(obj);
        end
    end
end

