function res = poissonNd_Neumann(f,lx)

%
%
%  res = poissonNd_Neumann(f,lx)
%
%

if(exist('lx','var')==0)
  lx                  = ones(1,5);
end


switch(numel(size(f)))
 case{2}
  res                 = poisson2d_Neumann(f,lx(1),lx(2));
 case{3}
  res                 = poisson3d_Neumann(f,lx(1),lx(2),lx(3));
 case{4}
  res                 = poisson4d_Neumann(f,lx(1),lx(2),lx(3),lx(4));
 otherwise
  error('Pb of dimension in poissonNd.m')
end