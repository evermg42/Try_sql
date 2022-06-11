function v = div_proj_1(u,rho0,rho1)

%
%
%   Copyright (c) 2012 Edouard Oudet
%

K  = @(X)pd_operator(X, +1);
KS = @(X)pd_operator(X, -1);
% 
% u = KS(u);
lx                      = ones(length(u.dim),1);

d                         = u.dim;
% u.M{1}(1:end-1:end,:,:)   = 0;
% u.M{2}(:,1:end-1:end,:)   = 0;
% u.M{3}(:,:,1)             = rho0;
% u.M{3}(:,:,end)           = rho1;
v                         = u;
du                        = div(u,lx);
  %% 3-D %%
p                       = poisson3d_Neumann(-du,lx(1),lx(2),lx(3));

v.M{1}(2:(end-1),:,:)   = v.M{1}(2:(end-1),:,:) - diff(p,[],1)*d(1)/lx(1);
v.M{2}(:,2:(end-1),:)   = v.M{2}(:,2:(end-1),:) - diff(p,[],2)*d(2)/lx(2);
v.M{3}(:,:,2:(end-1))   = v.M{3}(:,:,2:(end-1)) - diff(p,[],3)*d(3)/lx(3);
% v.M{1}(1:end-1:end,:,:)   = 0;
% v.M{2}(:,1:end-1:end,:)   = 0;
% v.M{3}(:,:,1)             = rho0;
% v.M{3}(:,:,end)           = rho1;

% v = K(v);
end