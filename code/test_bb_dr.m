%% Copyright (c) 2013 Gabriel Peyre, Nicolas papadakis, Edouard Oudet

%%
% Implements the minimization of the generalized BB energy using DR.
%     min_{U,V} J_epsilon(V) + i_I(U) + i_S(U,V)
% where

addpath('toolbox/');


N=32;P=N;Q=N;


d = [N,P,Q];
epsilon = 1e-6;
%exposant of the generalized cost
alpha= 1; % should be in [0;1];

% helpers
mynorm = @(a)norm(a(:));
mynorm_staggered = @(a)norm(a(:));
mymin = @(x)min(min(min(x(:,:,:,3))));
sum3 = @(a)sum(a(:));
[Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));
gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
normalize = @(u)u/sum(u(:));
obstacle=zeros(N,P,Q);
%%
% Load the data.

test = 'mixture';
test = 'obsession';
test = 'gaussian';
test = 'obstacle';
sigma =  .1;
rho = .000001; % minimum density value
switch test
    case 'mixture'
        f0 = normalize( rho + gaussian(.2,.3,sigma) );
        f1 = normalize( rho + gaussian(.6,.7,sigma*.7) + .6*gaussian(.7,.4,sigma*.7) );
        
    case 'gaussian'
        f0 = normalize( rho + gaussian(.2,.2,sigma) );
        f1 = normalize( rho + gaussian(.8,.8,sigma) );
        epsilon=min(f0(:));
    case 'obsession'
        sigma = 0.0354;
        rho=1e-6;
        f0 = normalize( rho + gaussian(.2,.2,sigma) );
        f1 = normalize( rho + gaussian(.8,.8,sigma) );
    case 'obstacle'
        Im=imread('Labyrinthe.png');
        N=size(Im,1);
        P=size(Im,2);
        Q=2*N;  
        [Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));
        gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
        d = [N,P,Q];
        obstacle=zeros(N,P,Q);
        sigma=0.04;
        f0 = gaussian(.08,.08,sigma);
        f1 = gaussian(.92,.92,sigma);
        
       
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

% prox operators
proxJeps    = @(U,V,gamma)deal(U, proxJ(V,gamma,epsilon,alpha,obstacle));
proxI       = @(U,V,gamma)deal(div_proj(U), V);
proxS       = @(U,V,gamma)interp_proj(U,V);
% functionals

J = @(w)sum3(  sum(w(:,:,:,1:2).^2,4) ./ (max((w(:,:,:,3)),max(epsilon,1e-10)).^alpha)  );
%%
% parallel-DR (PPXA) algorithm
Prox  = { proxJeps proxI proxS };

%Prox  = { proxJeps proxI proxS };
omega = [1 1 1]/3.;


% initialization
Xu = staggered(d); Xu.M{3} = f_init;
Xv = interp(Xu);
Yu_k = {Xu Xu Xu}; Yv_k = {Xv Xv Xv};
Zu_k=Yu_k;
Zv_k=Yv_k;

K = length(Prox);
% k=1: J
% k=2: div=0
% k=3: V=interp(U)

mu = 1.98; % should be in ]0,2[
gamma = 1./230.; % should be >0
niter = 1000;
Jlist = []; Constr = [];MinVal = [];


figure(1)
sel = round(linspace(1,Q+1,20));
tic
for i=1:niter
    progressbar(i,niter);
    for k=1:K
        % Zk = prox_{gamma*Fk/omegak}(Yk)
        [Zu_k{k},Zv_k{k}] = Prox{k}( Yu_k{k}, Yv_k{k}, gamma/omega(k) );
    end
    if 0
        mydisp = @(a)fprintf('Should be 0: %e\n', mynorm(a));
        mydisp( Zu_k{1}.M{3}(:,:,1) - f0 );
        mydisp( Zu_k{1}.M{3}(:,:,end) - f1 );
        % should be 0
        mydisp(div(Zu_k{2}));
        mydisp(Zv_k{3}-interp(Zu_k{3}));
    end
    % Z = sum_k omegak Zk
    Zu = staggered(d); Zv = zeros(N,P,Q,3);
    for k=1:K
        Zu = Zu + omega(k) * Zu_k{k};
        Zv = Zv + omega(k) * Zv_k{k};
    end
    for k=1:K
        % Yk = Y + mu*(2*Z-X-Zk)
        Yu_k{k} = Yu_k{k} + mu*(2*Zu - Xu - Zu_k{k});
        Yv_k{k} = Yv_k{k} + mu*(2*Zv - Xv - Zv_k{k});
    end
    % X = X + mu*(Z-X)
    Xu = Xu + mu*( Zu - Xu );
    Xv = Xv + mu*( Zv - Xv );
    % record energy
    Jlist(i)  = J(interp(div_proj(Xu)));
    Constr(i) = mynorm(div(Xu));
    MinVal(i)=mymin(interp(Xu));
    Xum=interp(Xu);
    
    
    
end
toc
%%
% Display the resulting density \(f(x,t)\) for \(t\) from 0 to 1.


clf;

close all
figure(1)

 if(max(obstacle(:)>0))
     obstacle(:,:,Q+1)=obstacle(:,:,Q); %increase the dimension for display
     U2=Xu.M{3};
     max_value=max(U2(:));
     I=obstacle>0;
     U2(I)= max_value;
     imageplot( mat2cell(U2(:,:,sel), N, P, ones(20,1)) , '', 2,3);axis equal
 else
     imageplot( mat2cell(Xu.M{3}(:,:,sel), N, P, ones(20,1)) , '', 2,3);
 end




figure(2)
clf;
subplot(2,1,1);
plot(Jlist);
title('J');
subplot(2,1,2);
plot(Constr); axis tight;
title('div=0 violation');
%disp(num2str(Constr(end)));



