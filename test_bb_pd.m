%%
% Implements the minimization of the generalized BB energy using PD.
%     min_{U} J_epsilon(K(U) + b) + i_I(U)
%  min_x F(K(U)) + G(U)   where   F(V)=J_epsilon(V + b)


%% Copyright (c) 2013 Gabriel Peyre, Nicolas papadakis, Edouard Oudet
close all
addpath('toolbox/');


N = 32; P = 32; Q = 32;

d = [N,P,Q];
epsilon = 1e-8;


%exposant of the generalized cost
alpha= 1; % should be in [0;1];


mynorm = @(a)norm(a(:));
dotp = @(a,b)sum(a(:).*b(:));
sum3 = @(a)sum(a(:));
[Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));
gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
normalize = @(u)u/sum(u(:));


%To represent the obstacles
obstacle=zeros(N,P,Q);

%%
% Load the data.
test = 'gaussian';
%test = 'obsession';
%test = 'mixture';

%test='obstacle';

sigma =  .1;
rho = 1e-12; % minimum density value
switch test
    case 'mixture'
        sigma = .1;
        sigma = .06;
        rho = .05; % minimum density value
        f0 = normalize( rho + gaussian(.2,.3,sigma) );
        f1 = normalize( rho + gaussian(.6,.7,sigma*.7) + .6*gaussian(.7,.4,sigma*.7) );
    case 'gaussian'
        f0 = normalize( rho + gaussian(.25,.25,sigma) );
        f1 = normalize( rho + gaussian(.75,.75,sigma) );
        epsilon=min(f0(:));
        
    case 'obsession'
        sigma = 0.0354;
        rho = 1e-4; % minimum density value
        % rho = .2;
        f0 = normalize( rho + gaussian(.25,.25,sigma) );
        f1 = normalize( rho + gaussian(.75,.75,sigma) );

     case 'obstacle'
%          sigma=0.04;
%              f0 = normalize( rho + gaussian(.2,.2,sigma) );
%         f1 = normalize( rho + gaussian(.8,.8,sigma) );
%         epsilon=min(f0(:));    
%          
%                   obstacle(15:24,13:15,:)=1;
%          obstacle(15:17,15:20,:)=1;
%          obstacle(20:20,25:25,:)=1;
        Im=imread('Labyrinthe.png');
        N=size(Im,1);
        P=size(Im,2);
        d = [N,P,Q];
        [Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));
        gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
        sigma=0.04;
        f0 = gaussian(.08,.08,sigma);
        f1 = gaussian(.92,.92,sigma);
        obstacle=zeros(N,P,Q);
       
        A=ones(N,P);
        I=Im(:,:,1)>0;
        A(I)=0;
        for i=1:Q,
            obstacle(:,:,i)=A;
        end

     
        I=obstacle(:,:,1)>0;
        f0(I)=0.;
        f1(I)=0.;
        f0=normalize(f0);
        f1=normalize(f1);
        
        epsilon=min(f0(:));
        
    otherwise
        error('Unknown');
end


clf; imageplot({f0 f1});

%%
% Initialization using linear interpolation of the densities.

t = repmat( reshape(linspace(0,1,Q+1), [1 1 Q+1]), [N P 1]);
f_init = (1-t) .* repmat(f0, [1 1 Q+1]) + t .* repmat(f1, [1 1 Q+1]);

% linear operator and adjoint
add_constraints = 1;

K  = @(X)pd_operator(X, +1);
KS = @(X)pd_operator(X, -1);
L = 1;
if 0
    % test for operator norm
    u = interp_adj(randn(N,P,Q,3));
    u = u*(1/norm(u)); e = [];
    for i=1:30
        v = KS(K(u));
        e(end+1) = dotp_stag(u,v);
        u = v * (1/norm(v));
    end
    L = e(end);
end

% to add into the J functional, J(K(U)+b) 
b = zeros(N,P,Q,3);
b(:,:,1,3) = f0/2; b(:,:,end,3) = f1/2;

% test for adjointness
U1 = interp_adj(randn(N,P,Q,3));
U2 = interp_adj(randn(N,P,Q,3));
if abs(dotp( K(U1), K(U2) ) - dotp_stag(KS(K(U1)), U2))>1e-5
    warning('Adjointness problem');
end

% prox operators
proxJeps  = @(V,gamma)proxJ(V,gamma,epsilon,alpha,obstacle);
proxF    = @(X,gamma)proxJeps(X+b,gamma)-b;
proxFS   = compute_dual_prox(proxF);
proxG    = @(U,gamma)pd_div_proj(U);
% functionals
J = @(V)sum3(  sum(V(:,:,:,1:2).^2,4) ./ (max((V(:,:,:,3)),max(epsilon,1e-10)).^alpha)  );  %pb interpU not >0  



%%
% Run the algorithm.

% initialization
U0 = staggered(d); U0.M{3} = f_init; 
% parameters
options.theta = 1.;
options.sigma=85;
options.tau = .99/(options.sigma*L);
options.niter = 2000;
mymin = @(x)min(min(min(x(:,:,:,3))));

options.report = @(U,V)struct( 'J', J(interp(U)), 'Constr', mynorm(div(U)), 'Min', mymin(interp(U)));%,...%


tic
[U,R,V] = perform_primal_dual(U0, K,  KS, proxFS, proxG, options);
toc
Jlist  = s2v(R,'J');
Constr = s2v(R,'Constr');
MinVal = s2v(R,'Min');


%%
% Display the resulting density \(f(x,t)\) for \(t\) from 0 to 1.

sel = round(linspace(1,Q+1,20));
V   = interp(U);
 figure;
 if(max(obstacle(:)>0))
     obstacle(:,:,Q+1)=obstacle(:,:,Q); %increase the dimension for display
     U2=U.M{3};
     max_value=max(U2(:));
     I=obstacle>0;
     U2(I)= max_value;
     imageplot( mat2cell(U2(:,:,sel), N, P, ones(20,1)) , '', 2,3);axis equal
 else
imageplot( mat2cell(U.M{3}(:,:,sel), N, P, ones(20,1)) , '', 2,3);axis equal
 end

 figure;
subplot(2,1,1);
plot(Jlist(10:end), '-'); axis tight;
title('J');
subplot(2,1,2);
plot((Constr(10:end)), '.-'); axis tight;
title('div=0 violation');




