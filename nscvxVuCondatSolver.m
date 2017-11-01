function [optsol, output] = nscvxVuCondatSolver(objFunc, linOper, x0, ...
                                                options, param, varargin)
% This implements Vu-Condat's algorithm to solve the following problem:
%
%          P^* := min_x { P(x) := g(Bx) + h(x) + r(x) },
%
% where g, h and r are convex, and r is Lipschitz gradient.
%
%
% Check input
if nargin < 5,       param   = [];                end
if nargin < 4,       options = PAPA_OptimSet([]); end
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
gxProxOper  = objFunc.gxProxOper;
hxProxOper  = objFunc.hxProxOper;
gxFunc      = objFunc.gxFunc;
hxFunc      = objFunc.hxFunc;
rxGradOper  = objFunc.rxGradOper; 
rxFunc      = objFunc.rxFunc;
Lipsr       = objFunc.rxLips;
Boper       = @(x) linOper.Boper( x, varargin{:});
BToper      = @(x) linOper.BToper(x, varargin{:});
LB_bar      = linOper.LB_bar;
gsProxOper  = @(x, gamma, varargin) ( x - gamma*gxProxOper(x/gamma, 1/gamma, varargin{:}) );
FxFunc      = @(x, varargin) gxFunc(Boper(x), varargin{:}) + hxFunc(x, varargin{:}) + rxFunc(x, varargin{:});        

% Set parameters.
if ~isempty(param)
    warning('PARAM must contain three values for tau, sigma and theta\n');
    tau_cur     = param.tau;
    sigma_cur   = param.sigma;
    theta_cur   = param.sigma;
else
    tau_cur     = 0.5/Lipsr;
    sigma_cur   = 0.5*Lipsr/LB_bar;
    theta_cur   = 1;
end

% Initialization.
y0          = zeros(size(Boper(x0)));
x_cur       = x0;
y_cur       = y0;

% The main loop.
for iter = 1:options.MaxIters
    
    % Evaluate the objective value.
    if options.isFxEval, fx_val = FxFunc(x_cur, varargin{:}); end
    
    % Start counting cummulative time.
    time_it = tic;
    
    % The main step.
    gradRx_cur   = rxGradOper(x_cur, varargin{:});
    x_tilde      = hxProxOper(x_cur - tau_cur * ( BToper(y_cur) + gradRx_cur ), tau_cur);
    y_tilde      = gsProxOper(y_cur + sigma_cur * Boper(2*x_tilde - x_cur), sigma_cur);
    x_next       = (1 - theta_cur)*x_cur + theta_cur*x_tilde;
    y_next       = (1 - theta_cur)*y_cur + theta_cur*y_tilde;
    
    % Compute the feasibility.
    abs_pfeas    = norm(y_next(:) - y_cur(:), 2);
    rel_pfeas    = abs_pfeas/max(norm(y_cur(:), 2), 1);
    
    % Compute the solution change.
    abs_schg     = norm(x_next(:) - x_cur(:), 2);
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
            y_cur        = y_next;
            break;
        end
    end
    
    % Move the the next iteration.
    x_cur        = x_next;
    y_cur        = y_next;
end
% End of the main loop.

% Finalization phase.
PAPA_finalization;

% Get the final solution.
optsol.x_opt  = x_cur;
optsol.y_opt  = y_cur;
optsol.fx_val = fx_val;

end

% PAPA v.1.0 by Quoc Tran-Dinh (quoctd@email.unc.edu)
% Copyright @ 2017 Department of Statistics and Operations Research (STOR)
%                The University of North Carolina at Chapel Hill (UNC)
% See the file LICENSE for full license information.