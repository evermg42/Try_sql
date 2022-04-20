function v = rigid_R(u)
% 
    v = zeros([u.dim 3]);
    [g1,g2,g3,g4] = gradient_x(u);
    constx = repmat([4;-4],[floor( (u.dim(1)-1)/2 ) u.dim(2:3)]);
    if mod(u.dim(1),2)==0
        constx = cat( 1, repmat(2,[1 u.dim(2:3)]), constx, repmat(4,[1 u.dim(2:3)]) );
    else
        constx = cat( 1, repmat(2,[1 u.dim(2:3)]), constx );
    end

    consty = repmat([4 -4],[u.dim(1) floor((u.dim(2)-1)/2) u.dim(3)]);
    if mod(u.dim(1),2)==0
        consty = cat( 2, repmat(2,[u.dim(1) 1 u.dim(3)]), consty, repmat(4,[u.dim(1) 1 u.dim(3)]) );
    else
        consty = cat( 2, repmat(2,[u.dim(1) 1 u.dim(3)]), consty );
    end
    
    %consty2 = 2 -2 2 -2 2 ...
    consty2 = cumprod([2;-ones(size(v,2),1)]);
    consty2 = repmat(consty2,[1 size(g2,2) size(g2,3)]);
    
    constx2 = cumprod([2 -ones(1,size(v,2))]);
    constx2 = repmat(constx2,[size(g3,1) 1 size(g3,3)]);
    
    for i = 1:u.dim(1)
        v(i,:,:,1) = sum(g1(i:-1:1,:,:).*constx(1:i,:,:),1);
    end
    
    for i = 1:u.dim(1)
        v(i,1:end-1,:,1) = v(i,1:end-1,:,1) + sum( g2(i:end,:,:).* consty2(i:end,:,:),1 );
        v(i,2:end,  :,1) = v(i,2:end,  :,1) - sum( g2(i:end,:,:).* consty2(i:end,:,:),1 );
    end


    for i = 1:u.dim(2)
        v(:,i,:,2) = sum(g4(:,i:-1:1,:).*consty(:,1:i,:),2);
    end

    for i = 1:u.dim(2)
        v(1:end-1,i,:,2) = v(1:end-1,i,:,2) + sum( g3(:,i:end,:).* constx2(:,i:end,:),2 );
        v(2:end,  i,:,2) = v(2:end,  i,:,2) - sum( g3(:,i:end,:).* constx2(:,i:end,:),2 );
    end

end