%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: DTxyz = divOperator(xy)
%%% PURPOSE:  Takes the vectorized and concatenated version of the
%%% gradients in both directions and then separates them and applies
%%% divergence operator.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DTxy = divOperator(Dxy, n)
    
    % Reshape the vector into an image.
    x           = Dxy(:,1:n);
    y           = Dxy(:,n+1:end);
    
    % Compute the divergence operator of x and y.
    DTxy        = [x(:,end) - x(:,1), -diff(x,1,2)];
    DTxy        = DTxy + [y(end,:) - y(1,:); -diff(y,1,1)];
end
