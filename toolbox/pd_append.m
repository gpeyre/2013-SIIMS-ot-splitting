function Y = pd_append(Y,u,v)

% aux function for PD algorithm

Y(:,:,1:2, end+1) = cat(3, u, v );