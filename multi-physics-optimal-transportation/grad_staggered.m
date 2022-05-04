function val = grad_staggered(u)
    d = u.dim;
    uxx = diff(u.M{1},1,1);
    uxy = diff(u.M{1},1,2);
    uyx = diff(u.M{2},1,1);
    uyy = diff(u.M{2},1,2);
    val = sum(uxx.^2,'all') + sum(uxy.^2,'all') + sum(uyx.^2,'all') + sum(uyy.^2,'all');
end