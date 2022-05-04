function [U,V] = interp_proj_v(U0,V0)
    U=U0;V=V0;
    d = U.dim;
    v = staggered(d(1:2));
    for ii=1:d(3)
        v.M{1} = U0.M{1}(:,:,ii);
        v.M{2} = U0.M{2}(:,:,ii);
        w = V0(:,:,ii,1:2);
        [v,w] = interp_proj(v,w);
        U.M{1}(:,:,ii) = v.M{1};
        U.M{2}(:,:,ii) = v.M{2};
        V(:,:,ii,1:2) = w;
    end
end