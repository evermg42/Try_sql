function u = grad_K(x,w)
    d = size(x);
    m1 = x(:,:,:,1);
    m2 = x(:,:,:,2);
    rho = x(:,:,:,3);
    v1 = w(:,:,:,1);
    v2 = w(:,:,:,2);
    u = zeros(size(x));
    u(:,:,:,1) = (m1 - rho .* v1);
    u(:,:,:,2) = (m2 - rho .* v2);
    u(:,:,:,3) = (v1.^2 + v2.^2) .* rho - (v1 .* m1 + v2 .* m2);
    