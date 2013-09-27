function V = interp(U)

%
%
% move from the staggered grid to the centered grid
%
%   V = interp(U);
%
%   U is staggered(dims)
%   V is of size [dims, d] where d=length(dims)
%
%   Copyright (c) 2012 Gabriel Peyre
%

switch(length(U.dim))
 case 2
   %% 2D field %%
   V          = cat(3, ...
                    U.M{1}(1:end-1,:) + U.M{1}(2:end,:), ...
                    U.M{2}(:,1:end-1) + U.M{2}(:,2:end) );
 case 3
   %% 3D field %%
   V          = cat(4, ...
                    U.M{1}(1:end-1,:,:) + U.M{1}(2:end,:,:), ...
                    U.M{2}(:,1:end-1,:) + U.M{2}(:,2:end,:), ...
                    U.M{3}(:,:,1:end-1) + U.M{3}(:,:,2:end));
 case 4
   %% 4D field %%
   V          = cat(5, ...
                    U.M{1}(1:end-1,:,:,:) + U.M{1}(2:end,:,:,:), ...
                    U.M{2}(:,1:end-1,:,:) + U.M{2}(:,2:end,:,:), ...
                    U.M{3}(:,:,1:end-1,:) + U.M{3}(:,:,2:end,:),...
		    U.M{4}(:,:,:,1:end-1) + U.M{4}(:,:,:,2:end));
 otherwise
  error('Undefined');
end
V             = V/2;
