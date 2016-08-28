function obj = projectonboundaryconditions(obj)

%
% CALL
%
% obj = projectonboundaryconditions(obj)
%
%



if(strcmp(obj.boundary.type,'Dirichlet'))	
  for k = 1:length(obj.M)
    obj.M{k}(obj.boundary.ind{k}) = obj.boundary.valb{k};
  end
end
