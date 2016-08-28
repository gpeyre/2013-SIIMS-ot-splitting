function [x,R] = perform_sdmm(x, ProxF, Mlist, MlistS, Qinv, options)

%
%
% perform_sdmm - simultaneous direction method of multipliers
%
%    [x,R] = perform_sdmm(x, {ProxF}, {Mlist}, {MlistS}, Qinv, options);
%
%   Solves min_x sum_i F_i(M_i*x)
%   where the F_i are convex proper functions with an easy to compute proximal operator.
%          
%
%   INPUTS  :
%   ProxF{i}(y,sigma) computes Prox_{sigma*F_i}(x)
%   Mlist{i}  the linear transforms
%   MlistS{i} the dual of the linear transforms (if isempty and Mlist{i} is numeric than a transpose matrix is used 
%   Qinv is the inverse of \sum_i Mlist{i}'*Mlist{i}. If isempty, a QR decomposition is used
%
%   options.niter is the number of iterations.
%   options.verb is for the diaplay of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS :
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   EX CALL : perform_sdmm(x, ProxF, Mlist, MlistS, Qinv, options)
%
%

m               = numel(ProxF);
niter           = getoptions(options, 'niter', 100);
verb            = getoptions(options, 'verb', 1);
gamma           = getoptions(options, 'gamma', 1);
report          = getoptions(options, 'report', @(x)0);

clear R;
countnum        = 0;
for k = 1:m
  if(isnumeric(Mlist{k}))
    Mfct{k}     = @(x) Mlist{k}*x;
    MfctS{k}    = @(x) Mlist{k}'*x;
    countnum    = countnum + 1;
  else
    Mfct{k}     = Mlist{k};
    MfctS{k}    = MlistS{k};
  end 
  y{k}          = Mfct{k}(x);
  z{k}          = y{k};
end
if(isnumeric(Qinv))
  if(and(isempty(Qinv),countnum==m)) % use a QR decomposition
    disp('Computation of the QR decomposition');
    Q           = Mlist{1}'*Mlist{1};
    for k = 2:m
      Q         = Q + Mlist{k}'*Mlist{k};
    end
    Rt          = qr(Q);
    Qinvfct     = @(x) Rt\(Rt'\(Q'*x(:)));
  else
    Qinvfct     = @(x) Qinv*x;
  end
else
  Qinvfct       = Qinv;
end


for i = 1:niter 
  xt            = 0;
  for k = 1:m
    xt          = xt +  MfctS{k}(y{k}-z{k});
  end
  x             = Qinvfct(xt);
  
  % record energies
  R(i)          = report(x);
  if(verb>1)
    fprintf('it %5d  %4.12f\n',i,R(i).J);
    %disp(num2str(x(:)'))
  elseif(verb)
    progressbar(i,niter);
  end
  for k = 1:m
    s{k}        = Mfct{k}(x);
    y{k}        = ProxF{k}(s{k}+z{k},gamma);
    z{k}        = z{k} + s{k} - y{k};
  end
end
