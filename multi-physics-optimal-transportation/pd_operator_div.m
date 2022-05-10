function Y = pd_operator_div(X, direction)

%
% pd_operator_div - linear operator for primal-dual scheme
sz = size(X);
Y = zeros(sz);
if(direction==1)
  % compute K
  for ii = 1:sz(3)
    Y(:,:,ii) = ( conv2(X(:,:,ii,1), [0 0 0;0  -sz(1) 0;0 sz(1) 0],'same' ) ) ...
              + ( conv2(X(:,:,ii,2), [0 0 0;0  -sz(2) sz(2);0 0 0],'same' ) );
  end
%   Y = zero_boundary(Y)
else
  % compute K^*
  for ii = 1:sz(3)
    Y(:,:,ii,1) = ( conv2(X(:,:,ii), [0 sz(1) 0;0 -sz(1) 0;0 0 0],'same' ) );
    Y(:,:,ii,2) = ( conv2(X(:,:,ii), [0 0 0;sz(2) -sz(2) 0;0 0 0],'same' ) );
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function X = zero_boundary(X)
%   
% 
% if(length(X.M)==2)
%   X.M{1}(1,:)     = 0; X.M{1}(end,:)     = 0;
%   X.M{2}(:,1)     = 0; X.M{2}(:,end)     = 0;
% elseif(length(X.M)==3)
%   X.M{1}(1,:,:)   = 0; X.M{1}(end,:,:)   = 0;
%   X.M{2}(:,1,:)   = 0; X.M{2}(:,end,:)   = 0;
%   X.M{3}(:,:,1)   = 0; X.M{3}(:,:,end)   = 0;
% elseif(length(X.M)==4)
%   X.M{1}(1,:,:,:) = 0; X.M{1}(end,:,:,:) = 0;
%   X.M{2}(:,1,:,:) = 0; X.M{2}(:,end,:,:) = 0;
%   X.M{3}(:,:,1,:) = 0; X.M{3}(:,:,end,:) = 0;
%   X.M{4}(:,:,:,1) = 0; X.M{4}(:,:,:,end) = 0;
% end
% end
end