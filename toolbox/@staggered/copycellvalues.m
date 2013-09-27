function C = copycellvalues(A,B)

%
% copy the values of the vector field B into the staggered class A (or vice versa) 
% (also the boundary values)
%

if(and(strcmp(class(A),'staggered'),strcmp(class(B),'cell')))
  C                             = A;
  for k = 1:length(C.dim)
    C.M{k}                      = B{k};
    C.boundary.valb{k}          = C.M{k}(C.boundary.ind{k});
  end
elseif(and(strcmp(class(A),'cell'),strcmp(class(B),'staggered')))
  C                             = B;
  for k = 1:length(C.dim)
    C.M{k}                      = A{k};
    C.boundary.valb{k}          = C.M{k}(C.boundary.ind{k});
  end
else
  warning('Strange classes in staggered.copycellvalues');
end
