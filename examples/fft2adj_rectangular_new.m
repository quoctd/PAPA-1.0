function y = fft2adj_rectangular_new(x, ind, m, n)
    y      = zeros([m, n]);
    y(ind) = x;
    y      = n*ifft2( ifftshift(y) );
end