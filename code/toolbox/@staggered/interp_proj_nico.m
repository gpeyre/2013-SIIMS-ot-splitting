function [U,V] = interp_proj_nico(U0,V0,b,A,A2,S2)


% projectect a staggered grid onto interpolation constraint
%
%   [U,V] = projection_staggered_interp(U0,V0);
%
%   U0 is staggered(dims) 
%   V0 is of size [dims, d] where d=length(dims)
%
% Copyright (c) 2012 Gabriel Peyre


% Compute projection of (u0,v0) on 
% C = {(u,v) in R^n x R^{n-1} \ v=S*u, u(1)=a=u0(1), u(2)=b=u0(end)}
% <=> A*[u;v] = g = [0...0;a;b]
% [u;v] = B*[u0;v0] + pA * g   where   B=Id-pA*A,  pA=A^+=A'*(AA')^{-1}


% projection matrices
persistent B;
persistent pA;
% dimension of V0, used to check wether need to update projection matrices
persistent dims; 

U = U0; V = V0;

d = length(U.dim);

for i=1:d
    
    %%% Extract and reshape %%%    
    switch d
        case 2
            u0 = U0.M{i};
            v0 = V0(:,:,i);    
           
            if i==2
                u0 = u0'; v0 = v0';
            end
        case 3
            u0 = U0.M{i};
            v0 = V0(:,:,:,i);
             b0=  b.M{i};
            if i==2
                u0 = permute(u0, [2 1 3]);
                v0 = permute(v0, [2 1 3]);
                b0 = permute(b0, [2 1 3]);
            elseif i==3
                u0 = permute(u0, [3 1 2]);
                v0 = permute(v0, [3 1 2]);
                b0 = permute(b0, [3 1 2]);
            end
            d_u0 = size(u0); d_v0 = size(v0); d_b0 = size(b0);
            u0 = reshape(u0, [d_u0(1) prod(d_u0(2:end))]);
            v0 = reshape(v0, [d_v0(1) prod(d_v0(2:end))]);
            b0 = reshape(b0, [d_b0(1) prod(d_b0(2:end))]);
        otherwise
            error('Only implemented for d=2 and d=3');
    end
    
    
    %%% Projection on the constraints %%%
    
   
    
    
    n = size(u0,1);
    k = size(u0,2);
    


    u=b0;
    C=S2'*v0-A2*b0;
    u(2:end-1,:)=A*(u0(2:end-1,:)+C(2:end-1,:));

    
    v=S2*u;
  
    
    %%% Reshape-back and insert %%%
    switch d
        case 2
            if i==2
                u = u'; v = v';
            end
            U.M{i} = u;
            V(:,:,i) = v;
        case 3
            u = reshape(u, d_u0);
            v = reshape(v, d_v0);
            if i==2
                u = permute(u, [2 1 3]);
                v = permute(v, [2 1 3]);
            elseif i==3
                u = permute(u, [2 3 1]);
                v = permute(v, [2 3 1]); 
            end
            U.M{i} = u;
            V(:,:,:,i) = v;
        otherwise
            error('Only implemented for d=2 and d=3');
            
    end
end

dims = size(V0);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% operator B of dimension (2*n-1)
function [pA,B] = projection_operator(n)

S = interp_operator(n);
A = [S, -eye(n); 1, zeros(1,2*n); zeros(1,n), 1, zeros(1,n)];
pA = pinv(A);
B = eye(2*n+1) - pA*A;

end

function S = interp_operator(N)

% interp_operator - load a midpoint interpolation operator 
%
%   S = interp_operator(N);
%
%   S as size (N,N+1), so it switches from size N+1 to size N
%
%   Copyrtight (c) 2012 Gabriel Peyre

[X,Y] = meshgrid(1:N+1,1:N+1);
S = eye(N+1) + ((X-Y)==1);
S = S(1:end-1,:)/2;

end