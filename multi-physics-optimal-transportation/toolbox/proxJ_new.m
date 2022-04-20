function V = proxJ_new(V0,gamm,epsilon,alpha,obstacle)
    %
    % proxJ - proximal operator of the J BB functional
    %
    %   V = proxJalpha(W,gamma,epsilon,alpha,c);
    %
    %   J(W) = sum_i c_i||m_i||^2/(f_i)^(alpha) + \Chi_{f_i>epsilon}
    %
    %%  W is assumed to be of dimension (N,P,Q,3)
    %
    %   Copyright (c) 2012 Gabriel Peyre
    %
    

    vs         = size(V0);
    d          = vs(end);
    %%[n,p,q,d]  = size(V0);
    %%if(d~=3)
    %%  error('Only works for d=3');
    %%end
    %%m0         = reshape(V0(:,:,:,1:2), [n*p*q 2]);
    %%f0         = reshape(V0(:,:,:,3  ), [n*p*q 1]);
    if(d==2)
      m0       = reshape(V0(:,:,1),       [prod(vs(1:(end-1))) 1]);
      f0       = reshape(V0(:,:,2),       [prod(vs(1:(end-1))) 1]);
    elseif(d==3)
      m0       = reshape(V0(:,:,:,1:2),   [prod(vs(1:(end-1))) 2]);
      f0       = reshape(V0(:,:,:,3),     [prod(vs(1:(end-1))) 1]);
    elseif(d==4)
      m0       = reshape(V0(:,:,:,:,1:3), [prod(vs(1:(end-1))) 3]);
      f0       = reshape(V0(:,:,:,:,4),   [prod(vs(1:(end-1))) 1]);
    else
      error('Only works for 2<= d <=4');
    end
    
    
    %Newton's method for finding the polynomial root
    
    % thisnorm = @(x) sum(x.^2,2);
    % thisnorm = @(x) cat(1, sum(x(1:end/2,:).^2*[4;1],2), sum(x(1+end/2:end,:).^2*[1;4],2));
    lambda = 4;
    mu = 0.01;
    A = [1 0;0 lambda];
    m0_new = V0(:,:,:,1:2);
    thisnorm = @(x) sum(x.^2,2);


    if alpha>0 && alpha<1 
        f=f0*0.+1.;  %Initialization f=1
    
    
        Numerateur = @(f)(f.^(2+alpha)+4.*gamm*f.^2-f0.*f.^(1+alpha)+4.*gamm^2*f.^(2-alpha)-4.*gamm*f0.*f-4.*gamm^2*f0.*f.^(1-alpha)-alpha*gamm*thisnorm(m0));
    
        Denominateur =@(f)(2+alpha)*f.^(1+alpha)+8.*gamm*f-(1+alpha)*f0.*f.^(alpha)+4.*(2-alpha)*gamm^2*f.^(1-alpha)-4.*gamm*f0-4.*(1-alpha)*gamm^2*f0./(f.^(alpha));
    
    
        Pnum=Numerateur(f);
        Pdiv=Denominateur(f);
    
        J=f>epsilon;
        k=0;
        reste=1.;
    
    
        %newton
        while k<50 && reste>1e-5   %threshold sufficient to have a correct estimation (taking thresholds<<1e-5 does not improve the results and slows down the process)   
        
       
    
            fm1=f;
            
            %the most different part proceed
            f(J)=(f(J)-(Pnum(J))./(Pdiv(J)));
              
    
            % not too small
            I=f<epsilon;
            f(I)       = epsilon;
            J=f>epsilon;
            
            Pnum=Numerateur(f);
            Pdiv=Denominateur(f);
    
            reste=max(abs(f(J)-fm1(J)));
     
            k=k+1;
        end
    
        
    elseif alpha==0
       f=f0; 
        
    else  %alpha==1
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!We only consider this
%         P          = [ones(length(f0),1), 4*gamm-f0, 4*gamm^2-4*gamm*f0, -gamm*thisnorm(m0) - 4*gamm^2*f0];
        P          = [ones(length(f0),1), 4*gamm-f0, 4*gamm^2-4*gamm*f0, -gamm*norm_new(m0_new,vs) - 4*gamm^2*f0];
        % roots
        R          = poly_root_new(P')';
        % positive root
        f          = real(R(:,1));
    end
    
    
    I          = f<epsilon;
    f(I)       = epsilon;
    I=obstacle>0;
    f(I)       = epsilon;
    % m          = m0./repmat(1+2*gamm./(f.^alpha), [1 (d-1)]);

    % tt = f.* m0;
    % ttt = (f .* eye(2) + gamm * (A + A'))';
    m = zeros(size(m0));
    for i = 1:size(m0,1)
        m(i,:) = (f(i).* m0(i,:)) / (f(i) .* eye(2) + gamm * (A + A'))';
    end
    V          = reshape([m f], vs); 

function X = norm_new(m0_new,vs)
    x_blocker = 1; y_blocker = 1;
    temp = x_blocker * m0_new(1:end/2,:,:,1).^2 + y_blocker * m0_new(1:end/2,:,:,2).^2;
    
    x_blocker = 1; y_blocker = 1;
    temp1 = x_blocker * m0_new(1+end/2:end,:,:,1).^2 + y_blocker * m0_new(1+end/2:end,:,:,2).^2;
    temp = cat(1, temp, temp1);
    X = reshape(temp, [prod(vs(1:(end-1))) 1]);
    
    