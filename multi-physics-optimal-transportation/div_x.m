function v = div_x(u)
% div - divergence operator
%   v = div(u);
%   if u is a 3-dimensional staggered grid, v is a 3-dimensional array.
% Copyright (c) 2022 Tao Ran

v = zeros(u.dim);
for k = 1:2
  v = v + diff(u.M{k},[],k);
end
