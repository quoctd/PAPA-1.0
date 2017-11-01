% FUNCTION: [l2nrm, err, iters] = PAPA_l2NormEval(nx, Aoper, AToper, maxiters, tol)
% PURPOSE: Compute the l2-norm of a matrix A^T*A by using Power method.
% 
% USAGE:
%    Inputs: 
%       nx          : The size of matrix A^T*A.
%       Aoper       : The linear operator A given by a function handle.
%       AToper      : The adjoint of A.
%       maxiters    : The maximum number of iterations (default, 30).
%       tol         : A given accuracy (default, 1e-7)
%    Outputs:
%       l2nrm       : The l2-norm of A^T*A.
%       err         : The output accuray
%       iters       : The number of iterations.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Department of Statistics and Operations Research
%    The University of North Carolina at Chapel Hill (UNC).
%    Date: 10.26.2017.
%    Last modified: 10.31.2017.
%    Contact: quoctd@email.unc.edu
%
function [l2nrm, err, iters] = PAPA_l2NormEval(nx, Aoper, AToper, maxiters, tol)

    if nargin < 5, tol = 1e-7; end
    if isempty(tol) || ~(tol > 0 && tol < 1), tol = 1e-7; end
    if nargin < 4, maxiters = 30; end
    if isempty(maxiters) || maxiters < 1, maxiters = 30; end

    % Initialization.
    if length(nx) == 1, nx = [nx, 1]; end
    z0     = ones(nx);
    q      = z0/norm(z0(:), 2); 
    relres = tol + 1; 
    z      = AToper(Aoper(q));

    % The main loop.
    for iter = 1:maxiters

        % Check the stopping criterion.
        if relres <= tol, break; end

        % Compute the new iteration.
        q       = z/norm( z(:), 2 ); 
        z       = AToper( Aoper(q) );
        lambda  = abs( sum(sum(q.*z)) );
        z2      = conj(z);
        q2      = z2/norm( z2(:), 2 ); 
        y1      = q2; 
        cosqy   = abs( sum(sum(y1.*q)) );

        if cosqy >= 5e-2
            relres = norm(z - lambda*q, 'fro')/cosqy;
        end
    end

    % Compute the norm of A'*A.
    l2nrm = abs(lambda);
    
    % Get the outputs.
    if nargout > 1, err   = relres; end
    if nargout > 2, iters = iter;   end
end

% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.