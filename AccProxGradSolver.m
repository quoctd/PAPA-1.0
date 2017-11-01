function [optsol, output] = AccProxGradSolver(objFunc, x0, options, varargin)
%
% This implement the accelerated proximal gradient algorithm to solve the
% following composite convex minimization problem
%
%                 F^* := min_x {F(x) := f(x) + g(x) },
%
% where f and g are convex, and f is Lipschitz gradient.
%

% Check input
if nargin < 3,       options = PAPA_OptimSet([]); end
if isempty(options), options = PAPA_OptimSet([]); end
options = PAPA_OptimSet(options);

% Get time started.
time1      = tic;

% Print the information and initialize the history list.
PAPA_printInfo();
PAPA_initHistory;

% Get problem size.
nx = objFunc.nx;
if length(nx) == 1, nx = [nx, 1]; end

% Get the problem input data.
fxGradOper  = objFunc.fxGradOper;
gxProxOper  = objFunc.gxProxOper;
fxFunc      = objFunc.fxFunc;
gxFunc      = objFunc.gxFunc;
fxLips      = objFunc.fxLips;

% Initialization.
x_cur       = x0;
x_hat       = x_cur;
t_cur       = 1;
stepSize    = 1.0/fxLips;

% The main loop.
for iter = 1:options.MaxIters
    
    % Evaluate the objective value.
    if options.isFxEval
        fx_val = fxFunc( x_cur, varargin{:} ) + gxFunc( x_cur, varargin{:} ); 
    end
    
    % Start counting cummulative time.
    time_it = tic;
        
    % The gradient step.
    gradFx_hat   = fxGradOper( x_hat, varargin{:} );
    x_next       = gxProxOper( x_hat - stepSize*gradFx_hat, stepSize, varargin{:} );
    
    % The accelerated step.
    t_next       = 0.5*( 1 + sqrt(1 + 4*t_cur.^2) );
    diff_x_cur   = x_next - x_cur;
    fact_cur     = (t_cur - 1)/t_next;
    t_cur        = t_next;
    x_hat        = x_next + fact_cur*diff_x_cur;
    
    % Compute the feasibility.
    abs_pfeas    = 0; rel_pfeas    = 0;
    
    % Compute the solution change.
    abs_schg     = norm(diff_x_cur(:), 2);
    rel_schg     = abs_schg/max(1, norm(x_cur(:), 2));
    
    % Print the iteration and save the iteration history.
    PAPA_printIteration;
    PAPA_saveHistory;
    
    % Check the stopping criterion.
    if options.isStoppingCond
        if rel_schg <= options.RelTolX && rel_pfeas <= options.RelTolFeas && iter > 1
            output.status = 'Convergence achieved';
            output.msg    = ['The serarch direction norm and the ', ...
                             'feasibility gap is below the desired threshold'];
            x_cur        = x_next;
            break;
        end
    end
    
    % Move the the next iteration.
    x_cur = x_next;
    
    % Perform Restarting if required.
    if options.isRestart
        if mod(iter, options.nRestart) == 0
            x_hat = x_cur;
            t_cur = 1;
        end
    end
end
% End of the main loop.

% Finalization phase.
PAPA_finalization;

% Get the final solution.
optsol.x_opt  = x_cur;
optsol.fx_val = fx_val;

end

% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.