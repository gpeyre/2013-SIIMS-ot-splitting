function C = mtimes(A,B)

%
% multiplication between staggered grid which respect boundary conditions
% or scalar multiplication
%

if(and(strcmp(class(A),'double'),strcmp(class(B),'staggered')))
  C                           = B;
  for k = 1:length(C.dim)
    C.M{k}                    = A*B.M{k};
  end
elseif(and(strcmp(class(A),'staggered'),strcmp(class(B),'double')))
  C                           = A;
  for k = 1:length(C.dim)
    C.M{k}                    = B*A.M{k};
  end
else
  C                           = A;
  for k = 1:length(C.dim)
    C.M{k}                    = A.M{k}.*B.M{k};
  end
end

% boundary conditions
C                             = projectonboundaryconditions(C);
