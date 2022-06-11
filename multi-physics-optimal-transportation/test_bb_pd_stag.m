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
obstacle=zeros(N,P,Q);
mysum = @(x)abs( sum(x(:))-1 );
%% Load the data
test = 'gaussian2';

sigma =  .1;
rho = 1e-12; % minimum density value
switch test
    case 'gaussian2'
        sigma = .07;
        f0 = normalize( rho + gaussian(.2,.4,sigma) + gaussian(.8,.6,sigma) );
        f1 = normalize( rho + gaussian(.4,.2,sigma) + gaussian(.6,.8,sigma) );
        epsilon=min(f0(:));
    otherwise
        error('Unknown');
end


% clf; imageplot({f0 f1});

%%
% Initialization using linear interpolation of the densities.

t = repmat( reshape(linspace(0,1,Q+1), [1 1 Q+1]), [N P 1]);
f_init = (1-t) .* repmat(f0, [1 1 Q+1]) + t .* repmat(f1, [1 1 Q+1]);

% linear operator and adjoint
% add_constraints = 1;

b = zeros(N,P,Q);
b(:,:,1) = f0; b(:,:,end) = f1;
bound = staggered(d);
K  = @(X)pd_operator_1(X, +1);
KS = @(X)pd_operator_1(X, -1);



L = 1;
% to add into the J functional, J(K(U)+b)

% test for adjointness
% U1 = interp_adj(randn(N,P,Q,3));
% U2 = interp_adj(randn(N,P,Q,3));
% if abs(dotp( K(U1), K(U2) ) - dotp_stag(KS(K(U1)), U2))>1e-5
%     warning('Adjointness problem');
% end

% prox operators
proxJeps  = @(V,gamma)proxJ_1(V,gamma,epsilon,alpha,obstacle);
proxG    = @(X,gamma)proxJeps(X,gamma);
proxF    = @(U,gamma)div_proj_1(U,f0,f1);
% proxFS   = @(y,sigma)y-sigma*proxF((1/sigma)*y,1/sigma);
proxFS   = compute_dual_prox(proxF);
% functionals
J = @(V)sum3(  sum(V(:,:,:,1:2).^2,4) ./ (max((V(:,:,:,3)),max(epsilon,1e-10)).^alpha)  );  %pb interpU not >0



%%
% Run the algorithm.

% initialization
y = staggered(d); y.M{3} = f_init;
U0 = y;
x1 = U0;
% parameters
theta = 1.;
sigma = 45;
tau = .02;% .99/(options.sigma*L);
mymin = @(x)min(min(min(x(:,:,:,3))));
% options.report = @(U,V)struct( 'J', J(interp(U)), 'Constr', mynorm(div(U)), 'Min', mymin(interp(U)));%,...%


tic

niter        = 20;
% y            = interp_adj(U0);
% y = zeros([d 3]);

%% iter£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡
for ii = 1:niter
    progressbar(ii,niter);
    xold       = U0;%
    y          = proxFS(y + sigma*x1,sigma);%

    U0         = KS(proxG(K(U0 - tau*y),tau));%

    x1         = U0 + theta*(U0-xold);%
    ysum = [sum(y.M{3}(:,:,1),'all') sum(y.M{3}(:,:,end),'all')]
    xsum = [sum(x1.M{3}(:,:,1),'all') sum(x1.M{3}(:,:,floor(end/2)),'all') sum(x1.M{3}(:,:,end),'all')]
    
    x = interp(U0);
    figure(2);
    subplot(2,2,1);quiver(x(:,:,1,2),x(:,:,1,1))
    subplot(2,2,2);quiver(x(:,:,floor(end/3),2),x(:,:,floor(end/3),1))
    subplot(2,2,3);quiver(x(:,:,floor(end*2/3),2),x(:,:,floor(end*2/3),1))
    subplot(2,2,4);quiver(x(:,:,end,2),x(:,:,end,1))
    figure(3);
    subplot(3,4,1);contour(x(:,:,1,3))
    subplot(3,4,2);contour(x(:,:,4,3));
    subplot(3,4,3);contour(x(:,:,7,3))
    subplot(3,4,4);contour(x(:,:,10,3))
    subplot(3,4,5);contour(x(:,:,13,3))
    subplot(3,4,6);contour(x(:,:,16,3))
    subplot(3,4,7);contour(x(:,:,19,3))
    subplot(3,4,8);contour(x(:,:,22,3))
    subplot(3,4,9);contour(x(:,:,25,3))
    subplot(3,4,10);contour(x(:,:,28,3))
    subplot(3,4,11);contour(x(:,:,30,3))
    subplot(3,4,12);contour(x(:,:,32,3))
end
toc




% Jlist  = s2v(R,'J');
% Constr = s2v(R,'Constr');
% MinVal = s2v(R,'Min');

%% Display the resulting density \(f(U0,t)\) for \(t\) from 0 to 1.

%  figure;
% subplot(2,1,1);
% semilogx(Jlist(10:end), '-'); axis tight;
% title('J');
% subplot(2,1,2);
% plot((Constr(10:end)), '.-'); axis tight;
% title('div=0 violation');
