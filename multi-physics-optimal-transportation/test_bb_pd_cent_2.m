clear;close all
addpath('toolbox/');
N = 32; P = 32; Q = 32;
d = [N,P,Q];
epsilon = 0;
%exposant of the generalized cost
alpha= 1; % should be in [0;1];
mynorm = @(a)norm(a(:));dotp = @(a,b)sum(a(:).*b(:));sum3 = @(a)sum(a(:));
[Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));
gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
normalize = @(u)u/sum(u(:));
obstacle=zeros(N,P,Q);
mysum = @(x)abs( sum(x(:))-1 );
%% Load the data
test = 'gaussian2';

sigma =  .1;
rho = 1e-12; % minimum density value
switch test
    case 'gaussian2'
        sigma = .07;
        f0 = normalize( rho + gaussian(.2,.5,sigma) + gaussian(.8,.5,sigma) );
        f1 = normalize( rho + gaussian(.4,.2,sigma) + gaussian(.6,.8,sigma) );
%         epsilon=min(f0(:));
    otherwise
        error('Unknown');
end

% clf; imageplot({f0 f1});

%%
K  = @(X)pd_operator_1(X, +1);
KS = @(X)pd_operator_1(X, -1);
% prox operators
proxJeps  = @(V,gamma)proxJ_1(V,gamma,epsilon,alpha,obstacle);
proxG    = @(X,gamma)proxJeps(X,gamma);
proxF    = @(U,gamma)div_proj_1(U);
% proxF    = @(U,gamma)div_proj_1(U,f0,f1);
% proxFS   = @(y,sigma)y-sigma*proxF((1/sigma)*y,1/sigma);
proxFS   = compute_dual_prox(proxF);
% functionals
J = @(V)sum3(  sum(V(:,:,:,1:2).^2,4) ./ (max((V(:,:,:,3)),max(epsilon,1e-10)).^alpha)  );  %pb interpU not >0



%%
% Run the algorithm.

% initialization
b = zeros(N,P,Q);
b(:,:,1) = f0; b(:,:,end) = f1;
U0 = zeros([d 3]);U0(:,:,:,3) = b;
U_ = U0;
y = KS(U_);
staggered(d);
% parameters
theta = 1.;sigma = 45;tau = .02;
mymin = @(x)min(min(min(x(:,:,:,3))));
% options.report = @(U,V)struct( 'J', J(interp(U)), 'Constr', mynorm(div(U)), 'Min', mymin(interp(U)));%,...%


tic

niter        = 500;
% y            = interp_adj(U0);
% y = zeros([d 3]);

rholist = [];
mlist = [];
Jlist = [];

%% iter£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡
for ii = 1:niter
    progressbar(ii,niter);
    xold       = U0;%
    y          = proxFS(y + sigma*KS(U_),sigma);%
    tp = K(y);
    tp(:,:,1,:) = 0;
    tp(:,:,end,:) = 0;
    x=U0;%U0old = sum(x(:,:,1,3),'all')
    U0         = proxG(U0 - tau*tp,tau);%

    U_         = U0 + theta*(U0-xold);%
%     ysum = [sum(y(:,:,1,3),'all') sum(y(:,:,end,3),'all')]
%     xsum = [sum(x1(:,:,1,3),'all') sum(x1(:,:,floor(end/2),3),'all') sum(x1(:,:,end,3),'all')]

    rholist = [rholist mynorm(U_(:,:,:,3)-xold(:,:,:,3))./mynorm(U_(:,:,:,3))]; 
    mlist = [mlist mynorm( U_(:,:,:,1:2)-xold(:,:,:,1:2))./mynorm(U_(:,:,:,1:2))];
    Jlist = [Jlist J(U_)];
end
toc

x = (U0);

figure(1);
subplot(3,4,1);contour(x(:,:,1,3));
subplot(3,4,2);contour(x(:,:,4,3));
subplot(3,4,3);contour(x(:,:,7,3));
subplot(3,4,4);contour(x(:,:,10,3));
subplot(3,4,5);contour(x(:,:,13,3));
subplot(3,4,6);contour(x(:,:,16,3));
subplot(3,4,7);contour(x(:,:,19,3));
subplot(3,4,8);contour(x(:,:,22,3));
subplot(3,4,9);contour(x(:,:,25,3));
subplot(3,4,10);contour(x(:,:,28,3));
subplot(3,4,11);contour(x(:,:,30,3));
subplot(3,4,12);contour(x(:,:,32,3));

figure(2);
subplot(2,2,1);quiver(x(:,:,1,2),x(:,:,1,1))
subplot(2,2,2);quiver(x(:,:,floor(end/3),2),x(:,:,floor(end/3),1))
subplot(2,2,3);quiver(x(:,:,floor(end*2/3),2),x(:,:,floor(end*2/3),1))
subplot(2,2,4);quiver(x(:,:,end,2),x(:,:,end,1))
    
figure(4);
loglog(rholist(10:end));hold on;
loglog(mlist(10:end));
legend('rho','m');

figure(5);
semilogx(Jlist(10:end));
legend('J');

% figure(6)
% sel = round(linspace(1,Q,20));
% imageplot( mat2cell(x(:,:,sel,3), N, P, ones(20,1)) , '', 2,3);