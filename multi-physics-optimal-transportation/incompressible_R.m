function v  = incompressible_R(u,alpha)
    
    % div_u = div(u);
    % mynorm = @(a)norm(a(:));
    % v = alpha * interp_adj(div_u);
    p = div(u);
    v = u;
    d = u.dim;

    v.M{1}(2:(end-1),:,:) = v.M{1}(2:(end-1),:,:) - alpha * diff(p,[],1)*d(1);
    v.M{2}(:,2:(end-1),:) = v.M{2}(:,2:(end-1),:) - alpha * diff(p,[],2)*d(2);
    v.M{3}(:,:,2:(end-1)) = v.M{3}(:,:,2:(end-1)) - alpha * diff(p,[],3)*d(3);
    
  
    
    lx = ones(length(u.dim),1);
    d                         = u.dim;
    du                        = div(u,lx);
    v                         = u;
    %% 3-D %%
    p                       = poisson3d_Neumann(-du,lx(1),lx(2),lx(3));
    v.M{1}(2:(end-1),:,:)   = v.M{1}(2:(end-1),:,:) - diff(p,[],1)*d(1)/lx(1);
    v.M{2}(:,2:(end-1),:)   = v.M{2}(:,2:(end-1),:) - diff(p,[],2)*d(2)/lx(2);
    v.M{3}(:,:,2:(end-1))   = v.M{3}(:,:,2:(end-1)) - diff(p,[],3)*d(3)/lx(3);