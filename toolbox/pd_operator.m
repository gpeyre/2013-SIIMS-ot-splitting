function Y = pd_operator(X, direction)

%
% pd_operator - linear operator for primal-dual scheme
%
%   Y = operator_pd(X, direction)
%
%   direction==+1 to compute K
%   direction==-1 to compute K^*
%
%   Copyright (c) 2012 Gabriel Peyre
%

if(direction==1)
  % compute K
  Y               = interp(zero_boundary(X));
else
  % compute K^*
  Y               = zero_boundary(interp_adj(X));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = zero_boundary(X)
  

if(length(X.M)==2)
  X.M{1}(1,:)     = 0; X.M{1}(end,:)     = 0;
  X.M{2}(:,1)     = 0; X.M{2}(:,end)     = 0;
elseif(length(X.M)==3)
  X.M{1}(1,:,:)   = 0; X.M{1}(end,:,:)   = 0;
  X.M{2}(:,1,:)   = 0; X.M{2}(:,end,:)   = 0;
  X.M{3}(:,:,1)   = 0; X.M{3}(:,:,end)   = 0;
elseif(length(X.M)==4)
  X.M{1}(1,:,:,:) = 0; X.M{1}(end,:,:,:) = 0;
  X.M{2}(:,1,:,:) = 0; X.M{2}(:,end,:,:) = 0;
  X.M{3}(:,:,1,:) = 0; X.M{3}(:,:,end,:) = 0;
  X.M{4}(:,:,:,1) = 0; X.M{4}(:,:,:,end) = 0;
end