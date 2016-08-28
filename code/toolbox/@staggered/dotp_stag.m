function d = dotp_stag(X,Y)

%
% dotp_stag - dot product between staggered arrays
%
%   d = dotp_stag(X,Y);
%
% Copyright (c) 2012 Gabriel Peyre
%

dp    = @(x,y) sum(x(:).*y(:));
d     = 0;
for i = 1:length(X.dim)
  d   = d + dp(X.M{i},Y.M{i});
end