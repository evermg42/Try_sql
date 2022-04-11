function v = div_proj_for_v(u,lx)

%
%
% div_proj - perform the projection on div=0 constraint
%
%   v = div_proj(u);
%
%   Copyright (c) 2012 Edouard Oudet
%

if(exist('lx','var')==0)
  lx                      = ones(length(u.dim),1);
end
d                         = u.dim;
du                        = div(u,lx);
v                         = u;
for ii = 1:d(3)
  p                       = poisson2d_Neumann(-du,lx(1),lx(2));
  v.M{1}(2:(end-1),:,ii)     = v.M{1}(2:(end-1),:,ii) - diff(p,[],1)*d(1)/lx(1);
  v.M{2}(:,2:(end-1),ii)     = v.M{2}(:,2:(end-1),ii) - diff(p,[],2)*d(2)/lx(2);
end