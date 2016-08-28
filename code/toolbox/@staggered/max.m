function C = max(A,B)

%
% C = max(A,B)
%
% maximum of a staggered grid. If nargin == 1, max is a scalar 
%
%

if(nargin==1)
  C                           = Inf;
  for k = 1:length(A.dim)
    C                         = max(A.M{k}(:),Inf);
  end
elseif(nargin==2)
  C                           = A;
  for k = 1:length(A.dim)
    C.M{k}                    = max(A.M{k},B.M{k});
  end
end


