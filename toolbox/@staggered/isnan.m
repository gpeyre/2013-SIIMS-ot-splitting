function [boolt,nbinf] = isnan(A)

%
% [boolt,nbnan] = isnanA)
%
% true if one NaN is found in A.M{:}
%
%


boolt                           = 0;
nbnan                           = 0;
for k = 1:length(A.dim)
  nbt                           = sum(isnan(A.M{k}(:)));
  if(nbt>0)
    nbnan                       = nbnan + nbt;
    boolt                       = 1;
  end
end
