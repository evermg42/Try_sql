function v  = incompressible_R(u)
%     v = u;
%     v.M{1} = cat(1, u.M{1}(2,:,:)-u.M{1}(1,:,:),diff(u.M{1},2,1),u.M{1}(end,:,:)-u.M{1}(end-1,:,:));
%     v.M{2} = cat(2, u.M{2}(:,2,:)-u.M{2}(:,1,:),diff(u.M{2},2,2),u.M{2}(:,end,:)-u.M{2}(:,end-1,:));
    d=u.dim;
    v=u;
    for ii = 1:v.dim(3)
        v1x = conv2(u.M{1}(:,:,ii),[-1;2;-1],'same');
        v2x = conv2(u.M{2}(:,:,ii),[-1 1;1 -1],'full');
        v.M{1}(:,:,ii) = v1x + v2x(:,2:end-1);
        v2y = conv2(u.M{2}(:,:,ii),[-1 2 -1],'same');
        v1y = conv2(u.M{1}(:,:,ii),[-1 1;1 -1],'full');
        v.M{2}(:,:,ii) = v2y + v1y(2:end-1,:);
    end
end