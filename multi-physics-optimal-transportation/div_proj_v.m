function v = div_proj_v(u)
    d = u.dim;
    v = staggered(d);
    temp = staggered(d(1:2));
    for tt = 1:d(3)
        temp.M{1} = u.M{1}(:,:,tt);
        temp.M{2} = u.M{2}(:,:,tt);
        temp = div_proj(temp);
        v.M{1}(:,:,tt) = temp.M{1};
        v.M{2}(:,:,tt) = temp.M{2};
    end
end