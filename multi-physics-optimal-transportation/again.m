%% Copyright (c) 2022 TAO Ran
% Implements the minimization
% min_{U,V} J_epsilon(V) + i_I(U) + i_S(U,V) + K(U,V)
clear;
addpath('toolbox/');

% helpers
Max = @(x)max(x,[],'all');
mynorm = @(a)norm(a(:));
mymin = @(x)min(min(min(x(:,:,:,3))));
sum3 = @(a)sum(a(:));
normalize = @(u)u/sum(u(:));

%% Load the initial data.
N=32;P=N;Q=N;
d = [N,P,Q];
[Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));
gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
obstacle=zeros(d);

sigma =  .05;
rho = .0; % minimum density value
f0 = normalize( rho + gaussian(.2,.6,sigma) + gaussian(.8,.4,sigma) + gaussian(.5,.5,sigma) );
f1 = normalize( rho + gaussian(.6,.2,sigma) + gaussian(.4,.8,sigma) + gaussian(.45,.6,sigma));


f0 = double(imread('shape_1.png'));
f1 = double(imread('shape_65.png'));
f0 = f0(:,:,1)/255;
f1 = f1(:,:,1)/255;
epsilon=0;
J = @(w)sum3(  sum(w(:,:,:,1:2).^2,4) ./ max( w(:,:,:,3),max(epsilon,1e-10) )  );

figure(1);
subplot(1,2,1);contour(f0);
subplot(1,2,2);contour(f1);
%% Initialization
t = repmat( reshape(linspace(0,1,Q+1), [1 1 Q+1]), [N P 1]);
finit = (1-t) .* repmat(f0, [1 1 Q+1]) + t .* repmat(f1, [1 1 Q+1]);
% finit(:,:,2,:) = zeros(size(finit(:,:,2,:)));

%% niter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
niter =2;
 lambda2 = 1000;
gamma = 1./230;
alphaJ= 1; % should be in [0;1];
alpha= 2000;
sel = round(linspace(1,Q+1,20));
new_sel = round(linspace(1,Q,20));

y0 = staggered(d);y0.M{3} = finit;
x0 = interp(y0);% the order: m1, m2, rho
q = staggered(d);
v0 = staggered(d);
w0 = interp(v0);

x=x0;y=y0;v=v0;w=w0;
old_x=x0;old_y=y0;old_v=v0;old_w=w0;

%part 1 prepare
sigma1 = 45;
tau1 = 0.02;
lambda = 0;
theta1 = .9;
%part 2 prepare
% sigma2 = ;
% proxCv    = @(U,V,gamma)deal(div_proj_v(U),V);
% proxS2     = @(U,V,gamma)interp_proj_v(U,V);
% proxG2    = {proxCv,proxS};

tt = x;
tt1= w;

L1 = 0; L2 = 0;

Jlist = [];
listw = [];
divxx = [];
gradv = [];
divxv = [];
%% Main

Type = "incompressible";
la = 1.98;
R = 0;
tic
for ii=1:niter
    progressbar(ii,niter);
    %% part 1
%     L1 = 2*lambda*max(1,max(w(:,:,:,1).^2 + w(:,:,:,2).^2,[],'all'))^2; %! F is L-Lipschitz continuous
%     for it=1:3
%         y = interp_adj(x0);
%         q = q + sigma1*y - sigma1 * div_proj((1/sigma1)*q + y);
%         old_x = x;
%         x = proxJ_K(x-tau1*interp(q),w,tau1,epsilon,lambda,obstacle);
%         x0 = x + theta1 * (x - old_x);
%         max(x0(:))
%     end
%     fprintf("(m1,m2,rho): D: %d %d %d norm: %d, %d, %d\n v should be %d %d\n", ...
%         mean(abs(x(:,:,:,1)-tt(:,:,:,1)),'all'), mean(abs(x(:,:,:,2)-tt(:,:,:,2)),'all'),mean(abs(x(:,:,:,3)-tt(:,:,:,3)),'all'), ...
%         mean(abs(x(:,:,:,1)),'all' ),                 mean(abs(x(:,:,:,2)),'all' ),                 mean(abs(x(:,:,:,3)),'all' ), ...
%         mean(abs(x(:,:,:,3)./x(:,:,:,1)),'all' ),     mean(abs(x(:,:,:,3)./x(:,:,:,2)),'all' ) );
%     tt = x;
%     %% part 2
%     
%     switch Type
%         case "div_0"
%             L2 = lambda2*max(x(:,:,:,3),[],'all').^2;
%             L2v = alpha*4*sqrt(2);
%         case "incompressible"
%             L2 = lambda2*max(x(:,:,:,3),[],'all').^2;
%             L2v = alpha*4*sqrt(2);
%         case "rigid"
%             L2 = lambda2*max(x(:,:,:,3),[],'all').^2;
%             L2v = alpha*4*sqrt(2);
%         case "translate"
%             L2 = lambda2*max(x(:,:,:,3),[],'all').^2;
%             L2v = alpha*4*sqrt(2);
%     end
% %     
%     tau = 0.6/(L2);
%     tauv = 0.6/(L2v);
%     zv = {v,v}; zw = {w,w};
%     delta_zv = zv;
%     delta_zw = zw;
%     for it=1:10
%         normw = 1;
%         cnt = 0;
%         while normw >= 1e-5 && cnt<=10
%             old_w =w;cnt = cnt+1;
%             % gradient of K
%             switch Type
%                 case "div_0"
%                     forward2 = lambda2*grad_K_w(x,w);
% %                     fprintf("%d %d\n",max(forward2(:)),tau);
%                     for typ = 1:2
%                         [delta_zv{typ},delta_zw{typ}] = proxG2{typ}( 2*v- zv{typ}, 2*w - zw{typ} - tau*forward2, gamma/omega(typ) );
%                         zv{typ} = zv{typ} + la * (delta_zv{typ}-v);
%                         zw{typ} = zw{typ} + la * (delta_zw{typ}-w);
%                     end
%                     v = .5 * zv{1}+.5 * zv{2};
%                     w = .5 * zw{1}+.5 * zw{2};
% 
%                 case "incompressible"
%                     R = incompressible_R(v);
%                     g = grad_K_w(x,w);
%                     forward2 = lambda2*g;
%                     forward2v = alpha*R;
%                     [v,w]  = proxS( v - tauv*forward2v, w - tau*forward2, gamma );
% %                     fprintf("g = %d R = %d div(v)=%d\n",mynorm(g(:)),mynorm([R.M{1}(:);R.M{2}(:)]),mynorm(div(v)));
%                 case "rigid"
%                     R = rigid_R(v);
%                     g = grad_K_w(x,w);
%                     forward2 = lambda2*g;
%                     forward2v = alpha*R;
%                     [v,w]  = proxS( v-tau*forward2v, w - tau*forward2, gamma );
%                 case "translate"
%                     R = translate_R(v);
%                     g = grad_K_w(x,w);
%                     forward2 = lambda2*g;
%                     forward2v = alpha*R;
%                     [v,w]  = proxS2( v - tauv*forward2v, w-tau*forward2 , gamma );
%             end
%             
%             
%             normw = mean(abs(w-old_w),'all');
%         end
% 
%     end
%     fprintf("(v1 v2 w1 w2): %d %d %d %d\n",mean(abs(v.M{1}),'all'),mean(abs(v.M{2}),'all'),mean(abs(w(:,:,:,1)),'all'),mean(abs(w(:,:,:,2)),'all') );
%     tttt = div_staggered(v);
%     fprintf("divxv = %d\n",tttt);
%     divxv = [divxv tttt];
%     tttt = grad_staggered(v);
%     fprintf("gradv = %d\n",tttt);
%     gradv = [gradv tttt];
%     tt1 = w;
end
toc
%%
% Display the result.
figure(2)
subplot(2,2,1)
quiver(w(:,:,1,2),w(:,:,1,1))
subplot(2,2,2)
quiver(w(:,:,floor(end/3),2),w(:,:,floor(end/3),1))
subplot(2,2,3)
quiver(w(:,:,floor(end*2/3),2),w(:,:,floor(end*2/3),1))
subplot(2,2,4)
quiver(w(:,:,end,2),w(:,:,end,1))

figure(3);
subplot(2,2,1)
quiver(x(:,:,1,2),x(:,:,1,1))
subplot(2,2,2)
quiver(x(:,:,floor(end/3),2),x(:,:,floor(end/3),1))
subplot(2,2,3)
quiver(x(:,:,floor(end*2/3),2),x(:,:,floor(end*2/3),1))
subplot(2,2,4)
quiver(x(:,:,end,2),x(:,:,end,1))
figure(7)
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
subplot(3,4,7)
contour(x(:,:,19,3))
subplot(3,4,8)
contour(x(:,:,22,3))
subplot(3,4,9)
contour(x(:,:,25,3))
subplot(3,4,10)
contour(x(:,:,28,3))
subplot(3,4,11)
contour(x(:,:,30,3))
subplot(3,4,12)
contour(x(:,:,32,3))