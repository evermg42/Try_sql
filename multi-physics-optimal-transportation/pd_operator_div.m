function Y = pd_operator_div(X, direction)
    % pd_operator_div - linear operator for primal-dual scheme
    sz = size(X);
    
    if(direction==1) % compute K
        Y = zeros(sz(1:end-1));
        for ii = 1:sz(3)
            Y(:,:,ii) = [conv2(X(:,:,ii,1), [sz(1);-sz(1)],'valid');zeros(1,sz(2))];
            Y(:,:,ii) = [conv2(X(:,:,ii,2), [sz(2) -sz(2)],'valid') zeros(sz(1),1)] + Y(:,:,ii);
        end

    else % compute K^*
        Y = zeros([sz 2]);
        for ii = 1:sz(3)
            Y(1:end-1,:,ii,1) = -X(1:end-1,:,ii);
            Y(2:end,  :,ii,1) =  X(1:end-1,:,ii)+Y(2:end,:,ii,1);
            Y(:,1:end-1,ii,2) = -X(:,1:end-1,ii);
            Y(:,2:end  ,ii,2) =  X(:,1:end-1,ii)+Y(:,2:end,ii,2);
        end
        Y(:,:,:,1) = sz(1)*Y(:,:,:,1);
        Y(:,:,:,2) = sz(2)*Y(:,:,:,2);
    end
end

