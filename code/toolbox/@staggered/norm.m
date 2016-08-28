function d = norm(X)

% norm - L2 norm of a staggered array
%
%   d = norm(X);
%
% Copyright (c) 2012 Gabriel Peyre

d = sqrt(dotp_stag(X,X));