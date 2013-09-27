function V = proxJ(V0,gamm,epsilon,alpha,obstacle)
%
% proxJ - proximal operator of the J BB functional
%
%   V = proxJalpha(W,gamma,epsilon,alpha,c);
%
%   J(W) = sum_i c_i||m_i||^2/(f_i)^(alpha) + \Chi_{f_i>epsilon}
%
%%  W is assumed to be of dimension (N,P,Q,3)
%
%   Copyright (c) 2012 Gabriel Peyre
%

vs         = size(V0);
d          = vs(end);
%%[n,p,q,d]  = size(V0);
%%if(d~=3)
%%  error('Only works for d=3');
%%end
%%m0         = reshape(V0(:,:,:,1:2), [n*p*q 2]);
%%f0         = reshape(V0(:,:,:,3  ), [n*p*q 1]);
if(d==2)
  m0       = reshape(V0(:,:,1),       [prod(vs(1:(end-1))) 1]);
  f0       = reshape(V0(:,:,2),       [prod(vs(1:(end-1))) 1]);
elseif(d==3)
  m0       = reshape(V0(:,:,:,1:2),   [prod(vs(1:(end-1))) 2]);
  f0       = reshape(V0(:,:,:,3),     [prod(vs(1:(end-1))) 1]);
elseif(d==4)
  m0       = reshape(V0(:,:,:,:,1:3), [prod(vs(1:(end-1))) 3]);
  f0       = reshape(V0(:,:,:,:,4),   [prod(vs(1:(end-1))) 1]);
else
  error('Only works for 2<= d <=4');
end


%Newton's method for finding the polynomial root

if alpha>0 && alpha<1  
    f=f0*0.+1.;  %Initialization f=1


    Numerateur = @(f)(f.^(2+alpha)+4.*gamm*f.^2-f0.*f.^(1+alpha)+4.*gamm^2*f.^(2-alpha)-4.*gamm*f0.*f-4.*gamm^2*f0.*f.^(1-alpha)-alpha*gamm*sum(m0.^2,2));

    Denominateur =@(f)(2+alpha)*f.^(1+alpha)+8.*gamm*f-(1+alpha)*f0.*f.^(alpha)+4.*(2-alpha)*gamm^2*f.^(1-alpha)-4.*gamm*f0-4.*(1-alpha)*gamm^2*f0./(f.^(alpha));


    Pnum=Numerateur(f);
    Pdiv=Denominateur(f);

    J=f>epsilon;
    k=0;
    reste=1.;


    %newton
    while k<50 && reste>1e-5   %threshold sufficient to have a correct estimation (taking thresholds<<1e-5 does not improve the results and slows down the process)   
    
   

        fm1=f;
  
        f(J)=(f(J)-(Pnum(J))./(Pdiv(J)));
          

        
        I=f<epsilon;
        f(I)       = epsilon;
        J=f>epsilon;
    Pnum=Numerateur(f);
    Pdiv=Denominateur(f);

        reste=max(abs(f(J)-fm1(J)));
 
        k=k+1;
    end

    
 elseif alpha==0
   f=f0; 
    
else  %alpha==1
    P          = [ones(length(f0),1), 4*gamm-f0, 4*gamm^2-4*gamm*f0, -gamm*sum(m0.^2,2) - 4*gamm^2*f0];
    % roots
    R          = poly_root_new(P')';
    % positive root
    f          = real(R(:,1));
    end


I          = f<epsilon;
f(I)       = epsilon;
I=obstacle>0;
f(I)       = epsilon;
m          = m0./repmat(1+2*gamm./(f.^alpha), [1 (d-1)]);

V          = reshape([m f], vs); 