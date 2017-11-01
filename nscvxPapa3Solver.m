% FUNCTION: [optsol, output] = nscvxPapaSolver(objFunc, linConstr, x0, y0, 
%                                              options, varargin)
% PURPOSE: An implementation of the PAPA algorithm for solving 
%          constrained convex optimization problems of the form:
%
%                    min_x g(x) + h(y) + r(y) s.t. -x + B*y in Kc,
%
%          where g, h and r are convex, and Kc is a given convex set.
%          Here, we do not assume g, h or r to be strongly convex.
%          But we assume that r is Lipschitz gradient.
%
% USAGE:
%    Inputs:
%      objFunc     : A structure that contains g and h, the proximal 
%                    operator of g, h and r, and additional fields if needed.
%      linConstr   : A structure to store the linear constraints with
%                    operators A and B, theirs adjoints and the projection 
%                    on the set Kc.
%      [x0, y0]    : An initial point
%      options     : Optional variables to control the algorithm.
%      varargin    : User defined variables if needed.
%
%    Outputs:
%      optsol      : A structure to contain the final solution.
%      output      : A structure to store other output information.
%
% REFFERENCES:
%    [1]. Q. Tran-Dinh, Proximal Alternating Pelnalized Algorithms for 
%         Fully Nonsmooth Constrained Convex Optimization, Manuscript, 2017.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Department of Statistics and Operations Research
%    The University of North Carolina at Chapel Hill (UNC).
%    Date: 10.14.2017.
%    Last modified: 10.14.2017.
%    Contact: quoctd@email.unc.edu
%
function [optsol, output] = nscvxPapa3Solver(objFunc, linConstr, x0, y0, ...
                                            options, varargin)
   
% Check the inputs.
if nargin < 2,           error('Not enough inputs');    end;
if nargin < 5,           options = PAPA_OptimSet([]);   end;
if nargin < 4,           y0      = [];                  end;
if nargin < 3,           x0      = [];                  end;
if isempty(options),     options = PAPA_OptimSet([]);   end;
if ~isstruct(linConstr), error('Invalid second input'); end
if ~isstruct(objFunc),   error('Invalid first input');  end

% Get time started.
time1       = tic;

% Print the information and initialize the history list.
PAPA_printInfo();
PAPA_initHistory;

% Get the size of the problem.
px          = objFunc.px;
py          = objFunc.py;
nc          = objFunc.nc;
if length(px) == 1, px = [px, 1]; end
if length(py) == 1, py = [py, 1]; end
if length(nc) == 1, nc = [nc, 1]; end

% Get the problem input data.
gxProxOper  = objFunc.gxProxOper;
hyProxOper  = objFunc.hyProxOper;
if isfield(objFunc, 'ryGradOper')
    ryGradOper = objFunc.ryGradOper; 
    ryFunc     = objFunc.ryFunc;
    Lipsr      = objFunc.ryLips;
end
gxFunc      = objFunc.gxFunc;
hyFunc      = objFunc.hyFunc;
Boper       = @(x) linConstr.Boper( x, varargin{:});
BToper      = @(x) linConstr.BToper(x, varargin{:});
dlStarProx  = linConstr.dlStarProx;
dlStarFunc  = linConstr.dlStarFunc;
beta1       = linConstr.beta1;

% Evaluate the l2-norm of B'*B.
LB_bar      = linConstr.LB_bar;
        
% Generate an initial point.
if isempty(x0), x0 = zeros(px); end
if isempty(y0), y0 = zeros(py); end

% Initialization phase.
x_cur       = x0;
y_cur       = y0;
y_hat       = y_cur;
lbd_dot     = zeros(nc);
tau         = 1;
beta        = beta1;
By_cur      = Boper(y_cur);
By_hat      = By_cur;

% The main loop.
for iter = 1:options.MaxIters
    
    % Evaluate the objective value if required.
    if options.isFxEval
        gx_val = gxFunc(x_cur, varargin{:}); 
        hy_val = hyFunc(y_cur, varargin{:});
        if isfield(objFunc, 'ryGradOper')
            ry_val = ryFunc(y_cur, varargin{:}); 
            fx_val = gx_val + hy_val + ry_val;
        else
            fx_val = gx_val + hy_val;
        end
        dl_val = dlStarFunc(-x_cur + By_cur, varargin{:});
    end
    
    % Start counting cumulative time.
    time_it = tic;
    
    % Compute the step-size tau.
    tau_next   = tau/(tau + 1);
    
    % Compute the step-size
    if isfield(objFunc, 'ryGradOper')
        betaL     = beta/(LB_bar + Lipsr*beta);
    else
        betaL      = beta/LB_bar;
    end
    
    % Compute the proximal operator of g and compute x_next.
    x_next     = gxProxOper( By_hat + beta*lbd_dot, beta, varargin{:} );
    
    % Compute the dual vector lbd_hat_star.
    lbd_hat_star = dlStarProx(lbd_dot + 1/beta * (-x_next + By_hat), 1/beta, varargin{:});
    
    % Compute the proximal operator of h and compute By_next.
    Btlbd_hats  = BToper(lbd_hat_star);
    if isfield(objFunc, 'ryGradOper')
        gradRy_hats = ryGradOper(y_hat, varargin{:});
        y_next  = hyProxOper( y_hat - betaL*(Btlbd_hats + gradRy_hats), betaL, varargin{:} );
    else
        y_next  = hyProxOper( y_hat - betaL*Btlbd_hats, betaL, varargin{:} );
    end
    By_next     = Boper(y_next);
    
    % Update vector x_hat and Ax_hat.
    diff_x_cur = x_next  - x_cur;
    diff_y_cur = y_next  - y_cur;
    y_hat      = y_next  + ((tau_next * (1 - tau))/tau) * diff_y_cur;
    By_hat     = By_next + ((tau_next * (1 - tau))/tau) * (By_next - By_cur);
    
    % Compute the solution change (absolute and relative changes).
    abs_schg   = norm([diff_x_cur(:); diff_y_cur(:)], 2);
    rel_schg   = abs_schg/max(1, norm([x_cur(:); y_cur(:)], 2));
    
    % Compute the feasibility gap (absotule and relative).
    abs_pfeas  = dl_val;
    if iter == 1, max_nrm_c1 = max(1, abs(abs_pfeas)); end
    rel_pfeas  = abs_pfeas/max_nrm_c1;
    
    % Print the iteration and save the iteration history.
    PAPA_saveHistory;
    PAPA_printIteration;
    
    % Check the terminated condition.
    if options.isStoppingCond,
        if rel_schg <= options.RelTolX && rel_pfeas <= options.RelTolFeas && iter > 1 
            output.status = 'Convergence achieved';
            output.msg    = ['The serarch direction norm and the ', ...
                             'feasibility gap is below the desired threshold'];
            x_cur         = x_next; y_cur = y_next;
            break;
        end
    end
    
    % Update the penalty parameter beta.
    beta  = (1 - tau_next)*beta;

    % Perform Restarting if required.
    if options.isRestart
        if mod(iter, options.nRestart) == 0
            y_hat     = y_next;
            By_hat    = By_next;
            lbd_dot   = lbd_hat_star;
            beta      = beta1;
            tau_next  = 1.0;
        end
    end
    
    % Assign variables to go to the next step.
    tau       = tau_next;
    x_cur     = x_next;
    y_cur     = y_next;
    By_cur    = By_next;
end
% End of the main loop.

% Finalization.
PAPA_finalization;

% Get the final solution.
optsol.x_opt   = x_cur;
optsol.y_opt   = y_cur;
optsol.lbd_opt = lbd_hat_star;
optsol.fx_val  = fx_val;

end

% PAPA v.1.0 by Quoc Tran-Dinh (quoctd@email.unc.edu)
% Copyright @ 2017 Department of Statistics and Operations Research (STOR)
%                The University of North Carolina at Chapel Hill (UNC)
% See the file LICENSE for full license information.