function b = LinearOperator(x, n, subrows, istrans)

    if istrans == 0
        Ab = dct(x);
        b  = Ab(subrows);
        return;
    end
    yfull          = zeros(n, 1);
    yfull(subrows) = x;
    b              = idct(yfull);
end