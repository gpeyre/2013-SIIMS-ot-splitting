function C = min(A,B)

%
% C = min(A,B)
%
% minimum of a staggered grid. If nargin == 1, min is a scalar 
%
%

if(nargin==1)
  C                           = Inf;
  for k = 1:length(A.dim)
    C                         = min(A.M{k}(:),Inf);
  end
elseif(nargin==2)
  C                           = A;
  for k = 1:length(A.dim)
    C.M{k}                    = min(A.M{k},B.M{k});
  end
end



