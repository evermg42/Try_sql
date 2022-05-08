function val = div_staggered(u)
    d = u.dim;
    uxx = diff(u.M{1},1,1);
    uyy = diff(u.M{2},1,2);
    val = sum((uxx+uyy).^2,'all');
end