function y = fft2fwd(x, ind, m, n)
    y = fftshift( fft2( x ) / m );
    y = y( ind );    
end
