function v = div_1(u,lx)

%
% div - divergence operator
%
%   v = div(u,lengths);
%
%   if u is a d-dimensional staggered grid, v is a d-dimensional array.
%
% Copyright (c) 2012 Edouard Oudet
%

if(exist('lx','var')==0)
  lx  = ones(length(u.dim),1);
end
d     = u.dim;
v     = 0;
for k = 1:length(d)
  v   = v + diff(u.M{k},[],k)*d(k)/lx(k);
end

end