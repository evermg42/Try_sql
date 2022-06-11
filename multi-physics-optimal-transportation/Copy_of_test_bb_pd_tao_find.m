%%
% Implements the minimization of the generalized BB energy using PD.
%     min_{U} J_epsilon(K(U) + b) + i_I(U)
%  min_x F(K(U)) + G(U)   where   F(V)=J_epsilon(V + b)
%% Copyright (c) 2013 Gabriel Peyre, Nicolas papadakis, Edouard Oudet
% close all;
clear;
addpath('toolbox/');
Constraint  = "divergence";
N = 32; P = 32; Q = 32;
d = [N,P,Q];
epsilon = 1e-8;
%exposant of the generalized cost

mynorm = @(a)norm(a(:));
dotp = @(a,b)sum(a(:).*b(:));
sum3 = @(a)sum(a(:));
[Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));
gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
normalize = @(u)u/sum(u(:));
%To represent the obstacles
obstacle=zeros(N,P,Q);

%% Load the data.
test = 'gaussian2';
rho = 0; % minimum density value
switch test
    case 'mixture'
        sigma = .06;
        rho = .05; % minimum density value
        f0 = normalize( rho + gaussian(.2,.3,sigma) );
        f1 = normalize( rho + gaussian(.6,.7,sigma*.7) + .6*gaussian(.7,.4,sigma*.7) );
    case 'gaussian'
        f0 = normalize( rho + gaussian(.25,.75,sigma) );
        f1 = normalize( rho + gaussian(.75,.25,sigma) );
        epsilon=min(f0(:));
    case 'gaussian2'
        sigma = .07;
        f0 = normalize( rho + gaussian(.2,.4,sigma) + gaussian(.8,.6,sigma) );
        f1 = normalize( rho + gaussian(.4,.2,sigma) + gaussian(.6,.8,sigma) );
        epsilon=min(f0(:));
     case 'shape'
        rho = 1e-4; % minimum density value
        f0 = double(imread('shape_1.png'));
        f0 = normalize(f0(:,:,1));
        f1 = double(imread('shape_65.png'));
        f1 = normalize(f1(:,:,1));
    otherwise
        error('Unknown');
end


% clf; imageplot({f0 f1});

%%
% Initialization using linear interpolation of the densities.

t = repmat( reshape(linspace(0,1,Q+1), [1 1 Q+1]), [N P 1]);
f_init = (1-t) .* repmat(f0, [1 1 Q+1]) + t .* repmat(f1, [1 1 Q+1]);

% linear operator and adjoint
add_constraints = 1;

K  = @(X)pd_operator(X, +1);
KS = @(X)pd_operator(X, -1);

switch  Constraint
    case "divergence"
        K2  = @(X)pd_operator_div(X, +1);
        KS2 = @(X)pd_operator_div(X, -1);
    case "rigidity"
        K2  = @(X)pd_operator_rigid(X, +1);
        KS2 = @(X)pd_operator_rigid(X, -1);
end
L = 1;

% to add into the J functional, J(K(U)+b)
b = zeros(N,P,Q,3);
b(:,:,1,3) = f0/2; b(:,:,end,3) = f1/2;

% test for adjointness
U1 = interp_adj(randn(N,P,Q,3));
U2 = interp_adj(randn(N,P,Q,3));
if abs(dotp( K(U1), K(U2) ) - dotp_stag(KS(K(U1)), U2))>1e-5
    warning('Adjointness problem');
end
U1 = randn([d 2]);
U2 = randn([d 2]);
sum( K2(U1).*K2(U2),'all' ) 
sum(KS2(K2(U1)).* U2,'all')
abs(sum( K2(U1).*K2(U2),'all' ) - sum(KS2(K2(U1)).* U2,'all'))
if abs(sum( K2(U1).*K2(U2),'all' ) - sum(KS2(K2(U1)).* U2,'all'))>1e-5
    warning('Adjointness problem');
end
% prox operators
w = zeros([d 2]);
% functionals

J = @(V)sum3(  sum(V(:,:,:,1:2).^2,4) ./ (max((V(:,:,:,3)),max(epsilon,1e-10)))  );  %pb interpU not >0
Kf = @(V,w).5*sum3(  (V(:,:,:,1:2)-V(:,:,:,3).*w(:,:,:,1:2)).^2  );
R = @(V) sum3(  K2(V).^2  );
% judgement
rholist = []; 
mlist = [];
qlist = [];

wlist = [];

Jlist = [];
Klist = [];
Rlist = [];
%% Run the algorithm.

% initialization
U0 = staggered(d); U0.M{3} = f_init;
U=U0;
U_ = U;
w_ = w;
% y = zeros([d 3]);
y = K(U0);
switch Constraint
    case "divergence"
        y2 = zeros(d);
    case "rigidity"
        y2 = zeros([d 4]);
end

% parameters

mu = 1e-4;
lambda = 10;

theta = 1.;
sigma=40;
tau = .02;

theta2 = 1.;
sigma2=4.5e-4;
tau2 = .01;

mymin = @(x)min(min(min(x(:,:,:,3))));


qold = y;q = y;
Vold = interp(U);V = interp(U);

proxG     = @(U,gamma)pd_div_proj(U);
proxF2   = @(X,gamma)X/(1+gamma*mu);
proxFS2  = compute_dual_prox(proxF2); 



Niter = 2000;
tic
for ii = 1:Niter
%     progressbar(ii,Niter);
    %% part 1
    %functions
    proxF     = @(X,gamma)proxJ_K(X+b,w,gamma,epsilon,lambda,obstacle)-b;
    proxFS    = compute_dual_prox(proxF);
    % update
    Uold       = U; 
    y          = proxFS(y + sigma*K(U_),sigma);
    U          = proxG(U - tau*KS(y),tau);
    U_         = U + theta*(U-Uold);

%     %% part 2
%     %functions
%     proxG2   = @(v,gamma)proxK(v,interp(U),gamma,lambda);
%     % update
%     wold       = w;
%     y2         = proxFS2(y2 + sigma2*K2(w_),sigma2);
%     w          = proxG2(w - tau2*KS2(y2),tau2);
%     w_         = w + theta2*(w-wold);
    %% check
    Vold = V;
    V = interp(U);
    rholist = [rholist mynorm(V(:,:,:,3)-Vold(:,:,:,3))./mynorm(V(:,:,:,3))]; 
    mlist = [mlist mynorm( V(:,:,:,1:2)-Vold(:,:,:,1:2))./mynorm(V(:,:,:,1:2))];
%     wlist = [wlist mynorm(w-wold)./mynorm(w)];
%     qold = q;q = y2;
%     qlist = [qlist mynorm(q-qold)./mynorm(q)];
    Jlist = [Jlist J(V)];
%     Klist = [Klist lambda*Kf(V,w)];
%     Rlist = [Rlist mu*R(w)];
end
toc
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

%  figure;
% subplot(2,1,1);
% plot(Jlist(10:end), '-'); axis tight;
% title('J');
% subplot(2,1,2);
% plot((Constr(10:end)), '.-'); axis tight;
% title('div=0 violation');

x = V;
figure;
subplot(3,4,1)
contour(x(:,:,1,3))
subplot(3,4,2)
contour(x(:,:,4,3))
subplot(3,4,3)
contour(x(:,:,7,3))
subplot(3,4,4)
contour(x(:,:,10,3))
subplot(3,4,5)
contour(x(:,:,13,3))
subplot(3,4,6)
contour(x(:,:,16,3))
% subplot(3,4,7)
% contour(x(:,:,19,3))
% subplot(3,4,8)
% contour(x(:,:,22,3))
% subplot(3,4,9)
% contour(x(:,:,25,3))
% subplot(3,4,10)
% contour(x(:,:,28,3))
% subplot(3,4,11)
% contour(x(:,:,30,3))
% subplot(3,4,12)
% contour(x(:,:,32,3))
% 
figure;
subplot(2,2,1)
quiver(w(:,:,1,2),w(:,:,1,1))
subplot(2,2,2)
quiver(w(:,:,floor(end/3),2),w(:,:,floor(end/3),1))
subplot(2,2,3)
quiver(w(:,:,floor(end*2/3),2),w(:,:,floor(end*2/3),1))
subplot(2,2,4)
quiver(w(:,:,end,2),w(:,:,end,1))

figure;
deltaw = KS2(y2);
subplot(2,2,1)
quiver(deltaw(:,:,1,2),deltaw(:,:,1,1))
subplot(2,2,2)
quiver(deltaw(:,:,floor(end/3),2),deltaw(:,:,floor(end/3),1))
subplot(2,2,3)
quiver(deltaw(:,:,floor(end*2/3),2),deltaw(:,:,floor(end*2/3),1))
subplot(2,2,4)
quiver(deltaw(:,:,end,2),deltaw(:,:,end,1))
figure;
loglog(rholist);hold on;
loglog(mlist);
loglog(qlist);
loglog(wlist);
legend('rho','m','q','w');
xlim([10,Niter]);
figure;
loglog(Jlist(10:end));hold on;
loglog(Klist);
loglog(Rlist);
xlim([10,Niter]);
legend('J','K','R');