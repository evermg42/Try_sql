function ProxFS = compute_dual_prox_1(ProxF)

% compute_dual_prox - compute the proximal operator of the dual

ProxFS = @(y,sigma)y-sigma*ProxF((1/sigma)*y,1/sigma);