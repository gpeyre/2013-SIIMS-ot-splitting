function C = minus(A,B)

%
% substraction between staggered grid which respect boundary conditions
%
%

C                             = A;
for k = 1:length(A.dim)
  C.M{k}                      = A.M{k} - B.M{k};
end

% boundary conditions
C                             = projectonboundaryconditions(C);
