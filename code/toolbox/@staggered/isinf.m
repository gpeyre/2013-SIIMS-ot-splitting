function [boolt,nbinf] = isinf(A)

%
% [boolt,nbinf] = isinfA)
%
% true if one Inf is found in A.M{:}
%
%


boolt                           = 0;
nbinf                           = 0;
for k = 1:length(A.dim)
  nbt                           = sum(isinf(A.M{k}(:)));
  if(nbt>0)
    nbinf                       = nbinf + nbt;
    boolt                       = 1;
  end
end
