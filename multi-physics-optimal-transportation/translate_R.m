function v = translate_R(u)%!
    v = zeros([u.dim 3]);
    [v1,v2]= gradient_x(u);
    v(:,1:end-1,:,1) = div_x(v1).^2;
    v(1:end-1,:,:,2) = div_x(v2).^2;
    