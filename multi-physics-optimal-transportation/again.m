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
%% Load the data.
N=20;P=N;Q=N;
d = [N,P,Q];
[Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));
gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
obstacle=zeros(d);

sigma =  .1;
rho = .000001; % minimum density value
f0 = normalize( rho + gaussian(.2,.5,sigma) + gaussian(.8,.5,sigma) );
f1 = normalize( rho + gaussian(.5,.8,sigma) + gaussian(.5,.2,sigma));
epsilon=min(f0(:));
J = @(w)sum3(  sum(w(:,:,:,1:2).^2,4) ./ max( w(:,:,:,3),max(epsilon,1e-10) )  );
%% Initialization
t = repmat( reshape(linspace(0,1,Q+1), [1 1 Q+1]), [N P 1]);
finit = (1-t) .* repmat(f0, [1 1 Q+1]) + t .* repmat(f1, [1 1 Q+1]);
finit(:,:,2,:) = zeros(size(finit(:,:,2,:)));

%% niter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
niter = 1000;
lambda = 1000; lambda2 = 1000;
gamma = 1./230;
alphaJ= 1; % should be in [0;1];
alpha= 2000;
sel = round(linspace(1,Q+1,20));
new_sel = round(linspace(1,Q,20));

y0 = staggered(d);y0.M{3} = finit;
x0 = interp(y0);% the order: m1, m2, rho

v0 = staggered(d);
w0 = interp(v0);

x=x0;y=y0;v=v0;w=w0;
old_x=x0;old_y=y0;old_v=v0;old_w=w0;

%part 1 prepare
proxJ1    = @(U,V,gamma)deal(div_proj(U),proxJ(V,gamma,epsilon,alphaJ,obstacle));
proxS     = @(U,V,gamma)interp_proj(U,V);
proxG1 = {proxJ1, proxS};
omega = [1/2 1/2];

%part 2 prepare
proxCv    = @(U,V,gamma)deal(div_proj(U),V);
proxG2    = {proxCv,proxS};

tt = x;
tt1= w;

L1 = 0; L2 = 0;

Jlist = [];
listw = [];
divxx = [];
%% Main

Type = "incompressible";
la = 1.98;
R = 0;
tic
for ii=1:niter
%     progressbar(ii,niter);
    %% part 1
    L1 = 2*lambda*max(1,max(w(:,:,:,1).^2 + w(:,:,:,2).^2,[],'all'))^2; %! F is L-Lipschitz continuous
    tau = .9/L1;
    

    zy = {y,y}; zx = {x,x};
    delta_zy = zy; delta_zx = zx;
    for it=1:10
        normx = 1;
        cnt =0;
        while normx>=1e-5 && cnt<=10
            cnt = cnt+1;
            old_x = x;
            forward = lambda * grad_K(x,w);
            for typ = 1:2
                [delta_zy{typ},delta_zx{typ}] = proxG1{typ}( 2*y - zy{typ}, 2*x - zx{typ} - tau*forward, gamma/omega(typ) );
            end
            for typ = 1:2
                zy{typ} = zy{typ} + la * (delta_zy{typ}-y);
                zx{typ} = zx{typ} + la * (delta_zx{typ}-x);
            end
            y =  omega(1) * zy{1}+ omega(2) * zy{2};
            x =  omega(1) * zx{1}+ omega(2) * zx{2};
            normx = max(abs(x-old_x),[],'all');
%             fprintf("%d\n",normx)
        end
    end
    Jlist(ii)  = J(interp(div_proj(y)));
%     figure(1);
%     imageplot( mat2cell(y.M{3}(:,:,1:Q+1), N, P, ones(Q+1,1)) , '', 2,3);
%     fprintf("LOL\n")
    fprintf("(m1,m2,rho): L: %d D: %d %d %d max: %d, %d, %d\n",L1,mean(abs(x(:,:,:,1)-tt(:,:,:,1)),'all'),mean(abs(x(:,:,:,2)-tt(:,:,:,2)),'all'),mean(abs(x(:,:,:,3)-tt(:,:,:,3)),'all'),max(x(:,:,:,1),[],'all'),max(x(:,:,:,2),[],'all'),max(x(:,:,:,3),[],'all'));
    tt = x;
    %% part 2
    
    switch Type
        case "div_0"
            L2 = lambda2*max(x(:,:,:,3),[],'all').^2;
        case "incompressible"
            L2 = lambda2*max(x(:,:,:,3),[],'all').^2 + alpha*4*sqrt(2);
        case "rigid"
            L2 = lambda2*max(x(:,:,:,3),[],'all').^2 + alpha*4*sqrt(2);
            
    end
    
    tau = 0.9/L2;
    zv = {v,v}; zw = {w,w};
    delta_zv = zv;
    delta_zw = zw;
    for it=1:10
        normw = 1;
        cnt = 0;
        while normw >= 1e-5 && cnt<=10
            old_w =w;cnt = cnt+1;
            % gradient of K
            switch Type
                case "div_0"
                    forward2 = lambda2*grad_K_w(x,w);
%                     fprintf("%d %d\n",max(forward2(:)),tau);
                    for typ = 1:2
                        [delta_zv{typ},delta_zw{typ}] = proxG2{typ}( 2*v- zv{typ}, 2*w - zw{typ} - tau*forward2, gamma/omega(typ) );
                        zv{typ} = zv{typ} + la * (delta_zv{typ}-v);
                        zw{typ} = zw{typ} + la * (delta_zw{typ}-w);
                    end
                    v = .5 * zv{1}+.5 * zv{2};
                    w = .5 * zw{1}+.5 * zw{2};

                case "incompressible"
                    R = incompressible_R(v);
                    g = grad_K_w(x,w);
                    forward2 = lambda2*g+alpha*R;
                    [v,w] = proxS( v, w - tau*forward2, gamma );
%                     fprintf("g = %d R = %d div_x(v)=%d\n",max(g(:)),max(R,[],'all'),max(div_x(v),[],'all'));
                    old_R=max(R,[],'all');
                case "rigid"
                    R = rigid_R(v);
                    g = grad_K_w(x,w);
                    forward2 = lambda2*g+alpha*R;
                    [v,w] = proxS( v, w - tau*forward2, gamma );
%                     fprintf("g = %d R = %d div_x(v)=%d\n",max(g(:)),max(R,[],'all'),max(div_x(v),[],'all'));
                    old_R=max(R,[],'all');
            end
        end
    end
%     [v11,v12,v21,v22] = gradient_x(v);
%     gradv =abs( [v11(:);v12(:);v21(:);v22(:)]);
%     fprintf("(w1,w2): L:%d D:%d %d max:%d %d grad_x(v)(sum,max)=(%d;%d)\n\n",L2,mean(abs(w(:,:,:,1)-tt1(:,:,:,1)),'all'),mean(abs(w(:,:,:,2)-tt1(:,:,:,2)),'all'),max(w(:,:,:,1),[],'all'),max(w(:,:,:,2),[],'all'),sum(gradv),max(gradv,[],'all'));
    fprintf("(w1,w2): L:%d D:%d %d max:%d %d div(v)=%d div_x(v)=%d\n\n",L2,mean(abs(w(:,:,:,1)-tt1(:,:,:,1)),'all'),mean(abs(w(:,:,:,2)-tt1(:,:,:,2)),'all'),max(w(:,:,:,1),[],'all'),max(w(:,:,:,2),[],'all'),mynorm(div(v)),mynorm(div_x(v)));
    divxx = [divxx mynorm(div(y))];
    tt1 = w;
end
toc
%%
% Display the result.
figure(2)
subplot(2,2,1)
quiver(w(end:-1:1,:,1,2),w(:,:,1,1))
subplot(2,2,2)
quiver(w(end:-1:1,:,floor(end/3),2),w(:,:,floor(end/3),1))
subplot(2,2,3)
quiver(w(end:-1:1,:,floor(end*2/3),2),w(:,:,floor(end*2/3),1))
subplot(2,2,4)
quiver(w(end:-1:1,:,end,2),w(:,:,end,1))

figure(3);
imageplot( mat2cell(y.M{3}(:,:,1:Q+1), N, P, ones(Q+1,1)) , '', 2,3);

figure(4)
subplot(2,2,1)
quiver(x(end:-1:1,:,1,2),x(:,:,1,1))
subplot(2,2,2)
quiver(x(end:-1:1,:,floor(end/3),2),x(:,:,floor(end/3),1))
subplot(2,2,3)
quiver(x(end:-1:1,:,floor(end*2/3),2),x(:,:,floor(end*2/3),1))
subplot(2,2,4)
quiver(x(end:-1:1,:,end,2),x(:,:,end,1))