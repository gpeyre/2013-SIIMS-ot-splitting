function U = interp_adj(V)

%
% interp_adj - adjoint of interp(@staggered)
%
%   U = interp_adj(V);
%
%   Satisfies <U,interpS(V)>=<interp(U),V>
%
% Copyright (c) 2012 Gabriel Peyre
%

d          = size(V); 
d          = d(1:end-1);
U          = staggered(d);
switch(length(d))
 case 2
   %% 2D field %%
   U.M{1}  = cat(1,V(1,:,1),V(1:end-1,:,1)+V(2:end,:,1),V(end,:,1));
   U.M{2}  = cat(2,V(:,1,2),V(:,1:end-1,2)+V(:,2:end,2),V(:,end,2));
 case 3
   %% 3D field %%
   U.M{1}  = cat(1,V(1,:,:,1),V(1:end-1,:,:,1)+V(2:end,:,:,1),V(end,:,:,1));
   U.M{2}  = cat(2,V(:,1,:,2),V(:,1:end-1,:,2)+V(:,2:end,:,2),V(:,end,:,2));
   U.M{3}  = cat(3,V(:,:,1,3),V(:,:,1:end-1,3)+V(:,:,2:end,3),V(:,:,end,3)); 
 case 4
   %% 4D field %%
   U.M{1}  = cat(1,V(1,:,:,:,1),V(1:end-1,:,:,:,1)+V(2:end,:,:,:,1),V(end,:,:,:,1));
   U.M{2}  = cat(2,V(:,1,:,:,2),V(:,1:end-1,:,:,2)+V(:,2:end,:,:,2),V(:,end,:,:,2));
   U.M{3}  = cat(3,V(:,:,1,:,3),V(:,:,1:end-1,:,3)+V(:,:,2:end,:,3),V(:,:,end,:,3)); 
   U.M{4}  = cat(4,V(:,:,:,1,4),V(:,:,:,1:end-1,4)+V(:,:,:,2:end,4),V(:,:,:,end,4)); 
 otherwise
  error('Only implemented for d=2, d=3 and d=4');
end
U          = U*.5;
