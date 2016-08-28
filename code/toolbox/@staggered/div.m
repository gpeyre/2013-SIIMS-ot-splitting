function v = div(u,lx)

%
% div - divergence operator
%
%   v = div(u,lengths);
%
%   if u is a d-dimensional staggered grid, v is a d-dimensional array.
%
% Copyright (c) 2012 Edouard Oudet
%

if(exist('lx','var')==0)
  lx  = ones(length(u.dim),1);
end
d     = u.dim;
v     = 0;
for k = 1:length(d)
  v   = v + diff(u.M{k},[],k)*d(k)/lx(k);
end












return
switch(length(d))
 case 2
  %% 2-D %%        
  v   = diff(u.M{1},[],1)*d(1) + diff(u.M{2},[],2)*d(2);
 case 3
  %% 3-D %%
  v   = diff(u.M{1},[],1)*d(1) + diff(u.M{2},[],2)*d(2) + diff(u.M{3},[],3)*d(3);
 otherwise 
  error('Only implemented for 2 <= d <= 4');
end