function res = poisson2d_Neumann(f,lx,ly)

%
%
%  res = poisson2d_Neumann(f,lx,ly)
%
%

if(exist('lx','var')==0)
  lx                  = 1;
end
if(exist('ly','var')==0)
  ly                  = 1;
end
N                     = size(f,1);
M                     = size(f,2);
hx                    = lx/N;
hy                    = ly/M;
dn                    = 0:1:N-1;
depn                  = 2*cos(pi*dn/N)-2;

dm                    = 0:1:M-1;
depm                  = 2*cos(pi*dm/M)-2;
denom2                = repmat(depn(:)/hx^2,1,M) + repmat(depm(:)'/hy^2,N,1);
denom2(denom2(:)==0)  = 1;
%for i=1:N
%   for j=1:M
%       denom(i,j) = depn(i)/hx^2 + depm(j)/hy^2;
%       if (denom(i,j)==0)
%           denom(i,j) = 1;
%       end
%   end
%end
%max(abs(denom(:)-denom2(:)))

fhat                  = mirt_dctn(f);
uhat                  = -(fhat)./denom2;
res                   = mirt_idctn(uhat);



