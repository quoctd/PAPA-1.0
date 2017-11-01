%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: Dxy = gradOperator(f_im)
%%% PURPOSE:  Define the gradient operator of a 2D image f. Vectorize and
%%% concatenate the gradients in both directions before giving the output.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dxy = gradOperator( f_im )
    
    % Compute the gradient of the images.
    Dx          = [diff(f_im, 1, 2), f_im(:, 1) - f_im(:, end)];
    Dy          = [diff(f_im, 1, 1); f_im(1, :) - f_im(end, :)];
    Dxy         = [Dx, Dy];
    
end
