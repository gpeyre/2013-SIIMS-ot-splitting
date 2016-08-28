function U = pd_div_proj(U)

% Just so some sanity check 

mynorm = @(x)norm(x(:));
mysum = @(x)abs( sum(x(:))-1 );

t = 1e-8;
if mynorm(U.M{1}(1,:,:))>t || mynorm(U.M{1}(end,:,:))>t || ...
   mynorm(U.M{2}(:,1,:))>t || mynorm(U.M{2}(:,end,:))>t || ...
   mysum( U.M{3}(:,:,1))>t || mysum( U.M{3}(:,:,end))>t 
    warning('Projection problem');
end

U = div_proj(U);