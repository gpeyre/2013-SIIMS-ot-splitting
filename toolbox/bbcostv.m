function Z = bbcostv(w,epsilon)

% 
%  N dimensiononal  generalization of 2D BB cost : sum3(sum(w(:,:,:,1:2).^2,4)./max(w(:,:,:,3),epsilon));  
%
%  Z = bbcostv(w) 
%

dim         = numel(size(w))-2;
sumv        = @(a) sum(a(:))/numel(a);
if(dim==1)
  Z         = sumv(w(:,:,1).^2./max(w(:,:,2),epsilon));  
elseif(dim==2)
  Z         = sumv(sum(w(:,:,:,1:2).^2,4)./max(w(:,:,:,3),epsilon));  
elseif(dim==3)
  Z         = sumv(sum(w(:,:,:,:,1:3).^2,5)./max(w(:,:,:,:,4),epsilon));  
end
