function V = proxJ_K(V0,w,tau,epsilon,lambda,obstacle)
%
% proxJ - proximal operator of the J BB functional
%
%   V = proxJalpha(W,gamma,epsilon,alpha,c);
%
%   J(W) = sum_i c_i||m_i||^2/(f_i)^(alpha) + \Chi_{f_i>epsilon}
%
%%  W is assumed to be of dimension (N,P,Q,3)
vs       = size(V0);
d        = vs(end);
m0       = reshape(V0(:,:,:,1:2),   [prod(vs(1:(end-1))) 2]);

m1       = reshape(V0(:,:,:,1),   [prod(vs(1:(end-1))) 1]);
m2       = reshape(V0(:,:,:,2),   [prod(vs(1:(end-1))) 1]);
f0       = reshape(V0(:,:,:,3),   [prod(vs(1:(end-1))) 1]);
w1       = reshape(w(:,:,:,1),    [prod(vs(1:(end-1))) 1]);
w2       = reshape(w(:,:,:,2),    [prod(vs(1:(end-1))) 1]);
absw2    = w1.^2 + w2.^2;
mTv      = m1.*w1+m2.*w2;

%Newton's method for finding the polynomial root

A = 1 + 4 * tau * lambda * ( 1 + .5*absw2 )  + 4 * tau^2 * lambda^2 * (1+absw2);
B = 2* tau - f0 + 4*tau^2*lambda - 4*tau*lambda*f0 -4*tau^2*lambda^2*f0 ...
    -2*tau*lambda*(mTv)-4*tau^2*lambda^2*(mTv) ...
    +4*tau^2*lambda*absw2+2*tau^3*lambda^2*absw2;
C = tau*( tau-2*f0-4*tau*lambda*(f0+mTv) + 2*tau^2*lambda*absw2 );
D = - 0.5*tau* ( (m1.^2 + m2.^2) + 2 * tau*f0 );
P = [A,B,C,D];
% roots
R          = poly_root_new(P')';
% positive root
f          = real(R(:,1));
I          = f<epsilon;
f(I)       = epsilon;
I=obstacle>0;
f(I)       = epsilon;
m          = (m0 + 2*lambda*tau * [f0 f0].*m0) ./ repmat(1+2*lambda*tau+tau./f, [1 (d-1)]);
V          = reshape([m f], vs);