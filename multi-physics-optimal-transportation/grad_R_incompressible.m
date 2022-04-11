function u = grad_R_incompressible(v)
    % alpha/lambda * grad_R
%     int_v = interp(v);
%     div_v = div(v);
    u = v;
%     sum_div_v = sum(div_v,3);
    
%     u.M{1}(2:end-1,:,:,:) = 2 * repmat(diff(sum_div_v,1,1),[1 1 32]).* diff(div_v,1,1);
%     u.M{2}(:,2:end-1,:,:) = 2 * repmat(diff(sum_div_v,1,2),[1 1 32]).* diff(div_v,1,2);
    d1 = diff(v.M{1}(:,:,:),1,1);
    d2 = diff(v.M{2}(:,:,:),1,2);
    u.M{1}(2:end-1,:,:) = (d1(2:end,:,:) + d2(1:end-1,:,:)).^2 - (d1(1:end-1,:,:) + d2(1:end-1,:,:)).^2;
    u.M{2}(:,2:end-1,:) = (d2(:,2:end,:) + d1(:,1:end-1,:)).^2 - (d2(:,1:end-1,:) + d1(:,1:end-1,:)).^2;