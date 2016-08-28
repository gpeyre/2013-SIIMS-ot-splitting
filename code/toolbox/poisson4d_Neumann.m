function res = poisson4d_Neumann(f,lx,ly,lz,lt)

%
%
%  res = poisson4d_Neumann(f,lx,ly,lz,lt)
%
%

if(exist('lx','var')==0)
  lx                  = 1;
end
if(exist('ly','var')==0)
  ly                  = 1;
end
if(exist('lz','var')==0)
  lz                  = 1;
end
if(exist('lt','var')==0)
  lt                  = 1;
end
N                     = size(f,1);
M                     = size(f,2);
R                     = size(f,3);
S                     = size(f,4);


hx                    = lx/N;
hy                    = ly/M;
hz                    = lz/R;
hw                    = lt/S;

dn                    = 0:1:N-1;
depn                  = 2*cos(pi*dn/N)-2;
dm                    = 0:1:M-1;
depm                  = 2*cos(pi*dm/M)-2;
dr                    = 0:1:R-1;
depr                  = 2*cos(pi*dr/R)-2;
ds                    = 0:1:S-1;
deps                  = 2*cos(pi*ds/S)-2;

denom2                = repmat(depn(:)/hx^2,[1,M,R,S]) + repmat(depm(:)'/hy^2,[N,1,R,S]) + ...
                        permute(repmat(depr(:)/hz^2,[1,N,M,S]),[2 3 1 4]) + ...
                        permute(repmat(deps(:)/hw^2,[1,N,M,R]),[2 3 4 1]);
denom2(denom2(:)==0)  = 1;
%for i = 1:N
%  for j = 1:M
%    for k = 1:R
%      for l = 1:S
%	denom(i,j,k,l)   = depn(i)/hx^2 + depm(j)/hy^2 + depr(k)/hz^2 + deps(l)/hw^2;
%	if (denom(i,j,k,l)==0)
%	  denom(i,j,k,l) = 1;
%	end
%      end
%    end
%  end
%end
%max(abs(denom(:)-denom2(:)))


fhat                  = mirt_dctn(f);
uhat                  = -(fhat)./denom2;
res                   = mirt_idctn(uhat);
