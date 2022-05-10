function v = proxK(v0,x,tau,lambda)
    m = x(:,:,:,1:2);
    rho = cat(4,x(:,:,:,3),x(:,:,:,3));
    v = (v0+tau*lambda*rho.*m)./(1+tau*lambda*rho.^2);
end