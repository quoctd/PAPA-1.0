function [optsol, output] = nscvxChambollePockSolver(objFunc, linOper, x0, ...
                                       options, varargin)
% This code implements Chambolle-Pock's algorithm to solve the following 
% composite convex optimization problem:
%
%                P^* := min_x { P(x) := f(x) + g(B*x) },
% 
% where f and g are convex and B is a linear operator.

% Get time started.
fprintf('******** The Chambolle-Pock primal-dual algorithm ********\n');
time1      = tic;

% Print the information and initialize the history list.
PAPA_printInfo();
PAPA_initHistory;

% Get data from inputs.
nx          = objFunc.nx;
proxFxOper  = objFunc.fxProxOper;
proxGxOper  = objFunc.gxProxOper;
proxFsOper  = @(x, gamma, varargin) ( x - gamma*proxFxOper(x/gamma, 1/gamma, varargin{:}) );
fxFunc      = objFunc.fxFunc;
gxFunc      = objFunc.gxFunc;
Aoper       = linOper.Aoper;
AToper      = linOper.AToper;
FxFunc      = @(x, varargin) fxFunc(Aoper(x), varargin{:}) + gxFunc(x, varargin{:});        
muhy        = 0;
if isfield(objFunc, 'muhy'), muhy = objFunc.muhy; end
if muhy > 0, fprintf('This is the strongly convex case\n'); end

% % Evaluate the l2-norm of AT*A.
LA_bar      = options.LA_bar;

% Set parameter.
tau_cur     = 0.5/sqrt(LA_bar);
sigma_cur   = tau_cur;
theta_cur   = 1;

% Initialization.
x_cur       = x0;
y0          = zeros(size(Aoper(x0)));
y_cur       = y0;
z_cur       = x0;

% The main loop.
for iter = 1:options.MaxIters
    
    % Start counting cummulative time.
    time_it = tic;
    
    % Evaluate the objective value.
    if options.isFxEval, fx_val = FxFunc(x_cur, varargin{:}); end
    
    % The main step.
    y_next       = proxFsOper(y_cur + sigma_cur * Aoper(z_cur), sigma_cur);
    x_next       = proxGxOper(x_cur  - tau_cur * AToper(y_next), tau_cur);
    theta_cur    = 1.0 / sqrt(1 + 2*muhy*tau_cur);
    z_cur        = x_next + theta_cur * (x_next - x_cur);

    % Update tau and sigma.
    tau_cur = theta_cur*tau_cur; sigma_cur = sigma_cur/theta_cur;
        
    % Compute the feasibility.
    abs_pfeas     = norm(y_next(:) - y_cur(:), 2);
    rel_pfeas     = abs_pfeas/max(norm(y_cur(:), 2), 1);
    
    % Compute the solution change.
    abs_schg      = norm(x_next(:) - x_cur(:), 2);
    rel_schg      = abs_schg/max(1, norm(x_cur(:), 2));
    
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

% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.
