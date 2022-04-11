%% Copyright (c) 2022 TAO Ran
% Implements the minimization
%     min_{U,V} J_epsilon(V) + i_I(U) + i_S(U,V) + K(U,V)
addpath('toolbox/');
% helpers
Max = @(x)max(x,[],'all');
mynorm = @(a)norm(a(:));
mymin = @(x)min(min(min(x(:,:,:,3))));
sum3 = @(a)sum(a(:));
normalize = @(u)u/sum(u(:));
J = @(w)sum3(  sum(w(:,:,:,1:2).^2,4) ./ max( w(:,:,:,3),max(epsilon,1e-10) )  );

%% Load the data.
N=7;P=N;Q=N;
d = [N,P,Q];
[Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));
gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
obstacle=zeros(d);

sigma =  .1;
rho = .000001; % minimum density value
f0 = normalize( rho + gaussian(.2,.2,sigma) );
f1 = normalize( rho + gaussian(.8,.8,sigma) );
epsilon=min(f0(:));

%% Initialization
t = repmat( reshape(linspace(0,1,Q+1), [1 1 Q+1]), [N P 1]);
finit = (1-t) .* repmat(f0, [1 1 Q+1]) + t .* repmat(f1, [1 1 Q+1]);


%% niter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
niter = 60;
lambda = 1000; lambda2 = 1;
gamma = 1./230;
alpha= 1; % should be in [0;1];
sel = round(linspace(1,Q+1,20));
new_sel = round(linspace(1,Q,20));

y0 = staggered(d);y0.M{3} = finit;
x0 = interp(y0);% the order: m1, m2, rho

v0 = staggered(d);
w0 = interp(v0);

x=x0;y=y0;v=v0;w=w0;
old_x=x0;old_y=y0;old_v=v0;old_w=w0;

%part 1 prepare
proxJ1    = @(U,V,gamma)deal(div_proj(U),proxJ(V,gamma,epsilon,alpha,obstacle));
proxS     = @(U,V,gamma)interp_proj(U,V);
proxG1 = {proxJ1, proxS};
omega = [1/2 1/2];

%part 2 prepare
proxCv    = @(U,V,gamma)deal(div_proj(U),V);
proxG2    = {proxCv,proxS};

tt = x;
tt1= w;
tt2 =y.M{3};

Jlist = [];
%% Main
tic
for ii=1:niter
%     progressbar(ii,niter);
    %% part 1
    la = 1.98;
    L = 2*lambda*max(1,max(w(:,:,:,1).^2 + w(:,:,:,2).^2,[],'all')); %! F is 1/be-Lipschitz continuous
    tau = .5/L;
    %%%
%     fprintf("ii:%d be:%d\n x:%d %d %d\n",ii,L,max(max(max(abs(x)))));
%     fprintf("w:%d %d %d\n",max(max(max(abs(w)))));
%     fprintf("y:%d\n",max(max(max(abs(y.M{3})))));
    fprintf("be: %d Dx: %d ",L,max(abs(x-tt),[],'all'));
    fprintf("Dw: %d\n",max(abs(w-tt1),[],'all'));
%     tt = x;
%     tt(:,:,:,3);
%     tt1= w;
%     tt2 =y.M{3};
    % 
    zy = {y,y}; zx = {x,x};
    delta_zy = zy; delta_zx = zx;
    for it=1:10
        normx = 1;
        while normx>=1e-4
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
%     imageplot( mat2cell(y.M{3}(:,:,1:Q+1), N, P, ones(Q+1,1)) , '', 2,3);
%     fprintf("LOL\n")

    %% part 2
    la = 1;
    L = lambda2*max(1,max(x(:,:,:,1).^2,[],'all')); %! F is 1/be-Lipschitz continuous
    tau = .5/L;
    zv = {v,v}; zw = {w,w};
    delta_zv = zv;
    delta_zw = zw;
    for it=1:10
        normw = 1;
        while normw >= 1e-4
            old_w =w;
            % gradient of K
            forward2(:,:,:,1) = (x(:,:,:,3)) .* w(:,:,:,1) - (x(:,:,:,3)).*x(:,:,:,1);
            forward2(:,:,:,2) = (x(:,:,:,3)) .* w(:,:,:,2) - (x(:,:,:,3)).*x(:,:,:,2);
            tpw(:,:,:,1) = w(:,:,:,1) - forward2(:,:,:,1);
            tpw(:,:,:,2) = w(:,:,:,2) - forward2(:,:,:,2);
            for typ = 1:2
                [delta_zv{typ},delta_zw{typ}] = proxG2{typ}( 2*v- zv{typ}, 2*w - zw{typ} - tau*forward2, gamma/omega(typ) );
            end
            for typ = 1:2
                zv{typ} = zv{typ} + la * (delta_zv{typ}-v);
                zw{typ} = zw{typ} + la * (delta_zw{typ}-w);
            end
            v = .5 * zv{1}+.5 * zv{2};
            w = .5 * zw{1}+.5 * zw{2};
            normw = max(abs(w-old_w),[],'all');
        end
    end

end
toc

%%
% Display the resulting density \(f(x,t)\) for \(t\) from 0 to 1.


figure(2);
imageplot( mat2cell(y.M{3}(:,:,1:Q+1), N, P, ones(Q+1,1)) , '', 2,3);
figure(3);
imageplot( {v.M{1}(:,1:P),v.M{2}(:,1:P)} );