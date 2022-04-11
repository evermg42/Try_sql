function [x, R, z] = GeneralizedForwardBackward1( gradF, proxGi, z, nIter, ga, la )

    n = length( proxGi );
    catDim = ndims( z );
    x = mean( z, catDim );
    N = numel( x );
    zi = zeros( size( x ) );

    for it=1:nIter
        forward = ga*gradF(x);
        for i=1:n
            idx = (i-1)*N+1:i*N;
            zi(:) = z(idx);	
            z(idx) = zi + la*( proxGi{i}( 2*x - zi - forward, n*ga ) - x );
        end
        x = mean( z, catDim );
    end

end %GeneralizedForwardBackward
