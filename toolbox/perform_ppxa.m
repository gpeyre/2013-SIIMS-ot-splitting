function [x,R] = perform_ppxa(x, ProxF, options)

%
%
% perform_ppxa - parallel proximal algorithm
%
%    [x,R] = perform_ppxa(x, {ProxF}, options);
%
%   Solves min_x sum_i F_i(x)
%   where the F_i are
%   convex proper functions with an easy to compute proximal operator.
%
%
%   INPUTS:
%   ProxF{i}(y,sigma) computes Prox_{sigma*F_i}(x)
%   options.niter is the number of iterations.
%   options.verb is for the diaplay of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%

m               = numel(ProxF);
mu              = getoptions(options, 'mu', 1);
niter           = getoptions(options, 'niter', 100);
verb            = getoptions(options, 'verb', 1);
gamma           = getoptions(options, 'gamma', 1);
omega           = getoptions(options, 'omega', ones(m,1)/m);
report          = getoptions(options, 'report', @(x)0);

clear R;
for k = 1:m
  y{k}          = x;
end
lambdat         = 1.3;
for i = 1:niter 
  % record energies
  R(i)          = report(x);
  if(verb>1)
    fprintf('it %5d  %4.12f\n',i,R(i).J);
    %disp(num2str(x(:)'))
  elseif(verb)
    progressbar(i,niter);
  end
  for k = 1:m
    pv{k}       = ProxF{k}(y{k},gamma/omega(k));
  end
  p             = omega(1)*pv{1};
  for k = 2:m
    p           = p + omega(k)*pv{k};
  end   
  for k = 1:m
    y{k}        = y{k} + lambdat*(2*p - x - pv{k});
  end
  x             = x + lambdat*(p-x);
end
