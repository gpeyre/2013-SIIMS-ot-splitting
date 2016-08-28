function v = div_proj(u,lx)

%
%
% div_proj - perform the projection on div=0 constraint
%
%   v = div_proj(u);
%
%   Copyright (c) 2012 Edouard Oudet
%

if(exist('lx','var')==0)
  lx                      = ones(length(u.dim),1);
end
d                         = u.dim;
du                        = div(u,lx);
v                         = u;
switch(length(d))
 case 2
  %% 2-D %%
  p                       = poisson2d_Neumann(-du,lx(1),lx(2));
  v.M{1}(2:(end-1),:)     = v.M{1}(2:(end-1),:) - diff(p,[],1)*d(1)/lx(1);
  v.M{2}(:,2:(end-1))     = v.M{2}(:,2:(end-1)) - diff(p,[],2)*d(2)/lx(2);
 case 3 
  %% 3-D %%
  p                       = poisson3d_Neumann(-du,lx(1),lx(2),lx(3));
  v.M{1}(2:(end-1),:,:)   = v.M{1}(2:(end-1),:,:) - diff(p,[],1)*d(1)/lx(1);
  v.M{2}(:,2:(end-1),:)   = v.M{2}(:,2:(end-1),:) - diff(p,[],2)*d(2)/lx(2);
  v.M{3}(:,:,2:(end-1))   = v.M{3}(:,:,2:(end-1)) - diff(p,[],3)*d(3)/lx(3);
 case 4
  %% 4-D %%
  p                       = poisson4d_Neumann(-du,lx(1),lx(2),lx(3),lx(4));
  v.M{1}(2:(end-1),:,:,:) = v.M{1}(2:(end-1),:,:,:) - diff(p,[],1)*d(1)/lx(1);;
  v.M{2}(:,2:(end-1),:,:) = v.M{2}(:,2:(end-1),:,:) - diff(p,[],2)*d(2)/lx(2);;
  v.M{3}(:,:,2:(end-1),:) = v.M{3}(:,:,2:(end-1),:) - diff(p,[],3)*d(3)/lx(3);;
  v.M{4}(:,:,:,2:(end-1)) = v.M{4}(:,:,:,2:(end-1)) - diff(p,[],4)*d(4)/lx(4);;
 otherwise
  error('Only implemented for 2 <= d <= 4');
end