function V = proxJ_1(V0,gamm,epsilon,alpha,obstacle)
%
% proxJ - proximal operator of the J BB functional
%
%   V = proxJalpha(W,gamma,epsilon,alpha,c);
%
%   J(W) = sum_i c_i||m_i||^2/(f_i)^(alpha) + \Chi_{f_i>epsilon}
%
%%  W is assumed to be of dimension (N,P,Q,3)
vs         = size(V0);
m0       = reshape(V0(:,:,:,1:2),   [prod(vs(1:(end-1))) 2]);
f0       = reshape(V0(:,:,:,3),     [prod(vs(1:(end-1))) 1]);
%alpha==1
P = [ones(length(f0),1), 4*gamm-f0, 4*gamm^2-4*gamm*f0,...
    -gamm*sum(m0.^2,2) - 4*gamm^2*f0];
% roots
R          = poly_root_new(P')';
% positive root
f          = real(R(:,1));


I          = f<epsilon;
f(I)       = epsilon;
I=obstacle>0;
f(I)       = epsilon;
m          = m0./repmat(1+2*gamm./(f.^alpha), [1 2]);

V          = reshape([m f], vs); 
mysum = @(x)abs( sum(x(:))-1 );
[mysum(V0(:,:,1,3)) mysum(V(:,:,1,3))];