function v  = incompressible_R(u)
% 
    v = zeros([u.dim 3]);
    divx = ( div_x(u) );
    
    constx = repmat([4;-4],[floor( (u.dim(1)-1)/2 ) u.dim(2:3)]);
    if mod(u.dim(1),2)==0
        constx = cat( 1, repmat(2,[1 u.dim(2:3)]), constx, repmat(4,[1 u.dim(2:3)]) );
    else
        constx = cat( 1, repmat(2,[1 u.dim(2:3)]), constx );
    end
    
    for i = 1:u.dim(1)
        v(i,:,:,1) = sum((divx(i:-1:1,:,:)).*constx(1:i,:,:),1);
    end
    
    consty = repmat([4 -4],[u.dim(1) floor((u.dim(2)-1)/2) u.dim(3)]);
    if mod(u.dim(1),2)==0
        consty = cat( 2, repmat(2,[u.dim(1) 1 u.dim(3)]), consty, repmat(4,[u.dim(1) 1 u.dim(3)]) );
    else
        consty = cat( 2, repmat(2,[u.dim(1) 1 u.dim(3)]), consty );
    end

    for i = 1:u.dim(2)
        v(:,i,:,2) = sum((divx(:,i:-1:1,:)).*consty(:,1:i,:),2);
    end
end
    