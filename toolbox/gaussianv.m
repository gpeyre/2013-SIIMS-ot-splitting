function Z = gaussianv(X,m,sigma)

%
%  exp(-((X{1}-m(1)).^2+(X{1}-m(2)).^2+...)/(2*sigma^2));
%
%

Z           = (X{1}-m(1)).^2;
for k = 2:numel(m)
  Z         = Z + (X{k}-m(k)).^2;
end
Z           = exp(-Z/(2*sigma^2));