function [L,e] = compute_operator_norm(A,n)

% compute_operator_norm - compute operator norm
%
%   [L,e] = compute_operator_norm(A,n);
%
%   Copyright (c) 2010 Gabriel Peyre

mynorm = @(x)norm(x(:));

if length(n)==1
    u = randn(n,1); u = u/mynorm(u);
else
    u = n;
    u = u/mynorm(u);
end
e = [];
for i=1:30
    v = A(u);
    e(end+1) = sum(u(:).*v(:));
    u = v/mynorm(v(:));
end
L = e(end);