function V = proxJ_new_1(V0,gamm,epsilon,alpha,obstacle)
    % proxJ - proximal operator of the J BB functional
    vs = size(V0);
    m0       = reshape(V0(:,:,:,1:2),   [prod(vs(1:(end-1))) 2]);
    f0       = reshape(V0(:,:,:,3),     [prod(vs(1:(end-1))) 1]);
    

    lambda = 4; mu = 0.1;
    A = zeros(2,2,size(V0,1),size(V0,2),size(V0,3));
    A(:,:,:, 1:end/2,     :) = repmat([lambda 0;0 1], [1 1 size(V0,1)  size(V0,2)/2 size(V0,3)]);
    A(:,:,:, 1+end/2:end, :) = repmat([mu 0;0 1], [1 1 size(V0,1)  size(V0,2)/2 size(V0,3)]);
    szA = size(A);
    A = reshape(A, [2 2 prod(szA(3:end))]);
    
    %assume alpha = 1
    f=f0*0.+1.;  %Initialization f=1

    
    A11=permute(A(1,1,:),[3 1 2]);
    A22=permute(A(2,2,:),[3 1 2]);
    m1=m0(:,1);
    m2=m0(:,2);
    % original function
    Numerateur = @(f)(f - f0) .* (2*gamm * A11 + f).^2 .* (2*gamm * A22 + f).^2 - ...
    gamm*A11.*(2 * gamm * A22 + f).^2 .* m1.^2 - gamm*A22.*(2 * gamm * A11 + f).^2 .* m2.^2;
    %its derivative
    Denominateur =@(f)5 * f.^4 + 12 * f.^2 * gamm .* (-A11 .* f0 - A22 .* f0 + A11.^2 * gamm + 4 * A11 .* A22 * gamm + ...
    A22.^2 * gamm) - 4 * f.^3 .* (f0 - 4 * (A11 + A22) * gamm) + ...
    4 * A11 .* A22 * gamm^2 .* (-4 * A22 .* f0 * gamm + 4 * A11 * gamm .* (-f0 + A22 * gamm) - m1.^2 - ...
    m2.^2) - 2 * f * gamm .* (4 * A11.^2 * gamm .* (f0 - 4 * A22 .* gamm) + ...
    A11 .* (16 * A22 .* f0 * gamm - 16 * A22.^2 * gamm^2 + m1.^2) + ...
    A22 .* (4 * A22 .* f0 * gamm + m2.^2));


    Pnum=Numerateur(f);
    Pdiv=Denominateur(f);

    J=f>epsilon;
    k=0;
    reste=1.;


    %newton
    
    while k<50 && reste>1e-5
        fm1=f;
        f(J)=(f(J)-(Pnum(J))./(Pdiv(J)));
        I=f<epsilon;
        f(I)       = epsilon;
        J=f>epsilon;
        Pnum=Numerateur(f);
        Pdiv=Denominateur(f);
        reste=max(abs(f(J)-fm1(J)));
        k=k+1;
    end
    I          = f<epsilon;
    f(I)       = epsilon;
    I=obstacle>0;
    f(I)       = epsilon;
    ttp = B(gamm, A, f);
    ttp = permute(ttp,[2 3 1]);
    m = zeros(size(m0));
    for ii = 1:size(ttp,3)
        m(ii,:) = f(ii) .* ( ttp(:,:,ii) \ m0(ii,:)' );
    end
    V          = reshape([m f], vs);

function X = B(gamm, A, f)
    tmp(1,:,:) = eye(2);
    tmp = repmat(tmp,[size(f,1) 1 1]);
    X = gamm * (permute(A, [3 2 1]) + permute(A, [3 1 2])) + f .* tmp;

function X = thisnorm(x, gamm, A, f)
    tmp = B(gamm, A, f);
    tmp = permute(tmp,[2 3 1]);
    X = zeros(32768,1);
    for ii = 1:size(tmp,3)
        t = tmp(:,:,ii) \ x(ii,:)';
        X(ii,:) = t' * A(:,:,ii) * t;
    end