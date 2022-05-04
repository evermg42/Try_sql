function v = translate_R(u)
    v = u;
    for ii = 1:v.dim(3)
    v.M{1}(:,:,ii) = conv2(u.M{1}(:,:,ii),[0 -1 0;-1 4 -1;0 -1 0],'same');
    v.M{2}(:,:,ii) = conv2(u.M{2}(:,:,ii),[0 -1 0;-1 4 -1;0 -1 0],'same');
%     v.M{2}(:,:,ii) = del2(u.M{2}(:,:,ii));
    end
end
    