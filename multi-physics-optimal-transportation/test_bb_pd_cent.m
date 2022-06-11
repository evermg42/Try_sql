clear;close all
addpath('toolbox/');


N = 32; P = 32; Q = 32;

d = [N,P,Q];
epsilon = 1e-8;
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
rho = 1e-12; % minimum density value
switch test
    case 'gaussian2'
        sigma = .07;
        f0 = normalize( rho + gaussian(.2,.5,sigma) + gaussian(.8,.5,sigma) );
        f1 = normalize( rho + gaussian(.4,.2,sigma) + gaussian(.6,.8,sigma) );
        epsilon=min(f0(:));
    otherwise
        error('Unknown');
end
Constraint = "divergence";

% clf; imageplot({f0 f1});

%%
switch  Constraint
    case "divergence"
        K2  = @(X)pd_operator_div(X, +1);
        KS2 = @(X)pd_operator_div(X, -1);
    case "rigidity"
        K2  = @(X)pd_operator_rigid(X, +1);
        KS2 = @(X)pd_operator_rigid(X, -1);
end
L = 1;


% U1 = randn([d 2]);
% U2 = randn([d 2]);
% sum( K2(U1).*K2(U2),'all' ) 
% sum(KS2(K2(U1)).* U2,'all')
% % abs(sum( K2(U1).*K2(U2),'all' ) - sum(KS2(K2(U1)).* U2,'all'))
% if abs(sum( K2(U1).*K2(U2),'all' ) - sum(KS2(K2(U1)).* U2,'all'))>1e-5
%     warning('Adjointness problem');
% end

% linear operator and adjoint
% add_constraints = 1;

K  = @(X)pd_operator_1(X, +1);
KS = @(X)pd_operator_1(X, -1);



L = 1;
% functionals
J = @(V)sum3(  sum(V(:,:,:,1:2).^2,4) ./ (max((V(:,:,:,3)),max(epsilon,1e-10)).^alpha)  );  %pb interpU not >0
Kf = @(V,w).5*sum3(  (V(:,:,:,1:2)-V(:,:,:,3).*w(:,:,:,1:2)).^2  );
R = @(V) sum3(  K2(V).^2  );


%%
% Run the algorithm.

% initialization
b = zeros(N,P,Q);
b(:,:,1) = f0; b(:,:,end) = f1;
u = zeros([d 3]);u(:,:,:,3) = b;
u_ = u;
y = KS(u);                            % ////////////////////////////

w = zeros([d 2]);
w_= w;
switch Constraint
    case "divergence"
        y2 = zeros(d);
    case "rigidity"
        y2 = zeros([d 4]);
end

% parameters
theta = 1.;sigma = 45;tau = .02;
mu = 5e-3;lambda = 30;
theta2 = 1.;sigma2=4.5e-4;tau2 = .01;

niter        = 300;

% prox operators
proxJeps  = @(V,gamma)proxJ_1(V,gamma,epsilon,alpha,obstacle);
proxF    = @(U,gamma)div_proj_1(U,f0,f1);
proxFS   = compute_dual_prox(proxF);

proxF2   = @(X,gamma)X/(1+gamma*mu);
proxFS2  = compute_dual_prox(proxF2); 



rholist = [];
mlist = [];
qold = K(y);q = qold;qlist = [];
wlist = [];
Jlist = [];
Klist = [];
Rlist = [];
tic
%% iter£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡
for ii = 1:niter
    progressbar(ii,niter);
    %% 1
    proxG     = @(X,gamma)proxJ_K(X,w,gamma,epsilon,lambda,obstacle);
    uold       = u;%
    y          = (proxFS((y + sigma*KS(u_)),sigma));%
    tp  =K(y);tp(:,:,1,:) = 0;tp(:,:,end,:) = 0;
    u         = proxG(u - tau*tp,tau);%
    u_         = u + theta*(u-uold);%

    rholist = [rholist mynorm(u_(:,:,:,3)-uold(:,:,:,3))./mynorm(u_(:,:,:,3))]; 
    mlist = [mlist mynorm( u_(:,:,:,1:2)-uold(:,:,:,1:2))./mynorm(u_(:,:,:,1:2))];
    Jlist = [Jlist J(u_)];
    %% 2
    proxG2   = @(v,gamma)proxK(v,u_,gamma,lambda);
    % update
    wold       = w;
    y2         = proxFS2(y2 + sigma2*K2(w_),sigma2);
    w          = proxG2(w - tau2*KS2(y2),tau2);
    w_         = w + theta2*(w-wold);

    wlist = [wlist mynorm(w-wold)./mynorm(w)];
    qold = q;q = tp;
    qlist = [qlist mynorm(q-qold)./mynorm(q)];
    
    %% judgement
    Klist = [Klist lambda*Kf(u_,w)];
    Rlist = [Rlist mu*R(w)];
end
toc


x = u_;
figure(1);
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

% figure(2);
% subplot(2,2,1);quiver(x(:,:,2,2),x(:,:,2,1))
% subplot(2,2,2);quiver(x(:,:,floor(end/3),2),x(:,:,floor(end/3),1))
% subplot(2,2,3);quiver(x(:,:,floor(end*2/3),2),x(:,:,floor(end*2/3),1))
% subplot(2,2,4);quiver(x(:,:,end-1,2),x(:,:,end-1,1))
% 
% figure(3);
% x = w;
% subplot(2,2,1);quiver(x(:,:,2,2),x(:,:,2,1))
% subplot(2,2,2);quiver(x(:,:,floor(end/3),2),x(:,:,floor(end/3),1))
% subplot(2,2,3);quiver(x(:,:,floor(end*2/3),2),x(:,:,floor(end*2/3),1))
% subplot(2,2,4);quiver(x(:,:,end-1,2),x(:,:,end-1,1))



figure(4);
loglog(rholist);hold on;
loglog(mlist);
loglog(qlist);
loglog(wlist);
legend('rho','m','q','w');



figure(5);
semilogx(Jlist);hold on;
semilogx(Klist);
semilogx(Rlist);
xlim([20,niter]);
legend('J','K','R');



% figure(6)
% 
% sel = round(linspace(1,Q,20));
% imageplot( mat2cell(x(:,:,sel,3), N, P, ones(20,1)) , '', 2,3);


