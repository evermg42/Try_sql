function v  = incompressible_R(u)
    v = u;
    v.M{1} = cat(1, u.M{1}(2,:,:)-u.M{1}(1,:,:),diff(u.M{1},2,1),u.M{1}(end,:,:)-u.M{1}(end-1,:,:));
    v.M{2} = cat(2, u.M{2}(:,2,:)-u.M{2}(:,1,:),diff(u.M{2},2,2),u.M{2}(:,end,:)-u.M{2}(:,end-1,:));
end