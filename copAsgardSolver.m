% FUNCTION: [optsol, output] = copAsgardSolver(objFunc, linConstr, x0, options, varargin)
% PURPOSE: An implementation of the ASGARD algorithm for solving 
%          constrained convex optimization problems of the form:
%                    min_x f(x) s.t. A*x = c.
%
% USAGE:
%    Inputs:
%      objFunc     : A structure that contains f, the proximal operator of f
%                    and additional fields if needed.
%      linConstr   : A structure to store the linear constraints with
%                    operators A, its adjoints and vector c.
%      x0          : An initial point
%      options     : Optional variables to control the algorithm.
%      varargin    : User defined variables if needed.
%
%    Outputs:
%      optsol      : A structure to contain the final solution.
%      output      : A structure to store other output information.
%
% REFFERENCES:
%    [1]. Q. Tran-Dinh, O. Fercoq, and V. Cevher, A Smoothing Primal-Dual 
%         Optimization Framework for Fully Nonsmooth Constrained Convex 
%         Optimization, Manuscript, 2016.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Department of Statistics and Operations Research
%    The University of North Carolina at Chapel Hill (UNC).
%    Joint work with Olivier Fercoq, and Volkan Cevher.
%    Date: 08.26.2016.
%    Last modified: 08.06.2016.
%    Contact: quoctd@email.unc.edu
%
function [optsol, output] = copAsgardSolver(objFunc, linConstr, x0, options, varargin)
   
% Check the inputs.
if nargin < 2,        error('Not enough inputs');    end;
if nargin < 4,        options = PAPA_OptimSet([]); end;
if nargin < 3,        x0      = [];                  end;
if isempty(options),  options = PAPA_OptimSet([]); end;
if ~isstruct(linConstr), error('Invalid second input');  end
if ~isstruct(objFunc),   error('Invalid first input');   end

% Get time started.
time1       = tic;
fprintf('******* ASGARD Solver *******\n');

% Print the information and initialize the history list.
PAPA_printInfo();
PAPA_initHistory;

% Get data from inputs.
nx          = objFunc.nx;
fxProxOper  = objFunc.fxProxOper;
fxFunc      = objFunc.fxFunc;
Aoper       = @(x) linConstr.Aoper( x, varargin{:});
AToper      = @(x) linConstr.AToper(x, varargin{:});
cb          = linConstr.cb;
mx          = length(Aoper(x0));
max_nrm_c1  = max(norm(cb, 2), 1);
gyStarProx  = linConstr.gyStarProx;
gyStarFunc  = linConstr.gyStarFunc;
beta1       = options.beta1;
fprintf(' ******* One-Loop Restarting ASGARD ******* \n');

% Evaluate the l2-norm of AT*A.
LA_bar     = options.LA_bar;
        
% Generate an initial point.
if isempty(x0), x0 = zeros(nx, 1); end

% Initialization phase.
x_cur       = x0;
x_hat       = x_cur;
y_dot       = zeros(mx, 1);
tau         = 1;
beta        = beta1;
Ax_cur      = Aoper(x_cur);
Ax_hat      = Ax_cur;

% The main loop.
for iter = 1:options.MaxIters
    
    % Start counting cumulative time.
    time_it = tic;
    
    % Evaluate the objective value if required.
    if options.isFxEval
       fx_val = fxFunc(x_cur, varargin{:}); 
       gy_val = gyStarFunc(Ax_cur, varargin{:});
    end
    
    % Compute the step-size tau.
    tau_next   = tau/(tau + 1);
    
    % Compute the dual vector y_hat_star.
    y_hat_star = gyStarProx(y_dot + 1/beta * Ax_hat, 1/beta);
    
    % Compute the proximal operator of f and compute Ax_next.
    betaL      = beta/LA_bar;
    Aty_hats   = AToper(y_hat_star);
    
    x_next     = fxProxOper( x_hat - betaL*Aty_hats, betaL, varargin{:} );
    Ax_next    = Aoper(x_next);
    
    % Update vector x_hat and Ax_hat.
    x_hat      = x_next  + ((tau_next * (1 - tau))/tau) * (x_next  - x_cur);
    Ax_hat     = Ax_next + ((tau_next * (1 - tau))/tau) * (Ax_next - Ax_cur);
    
    % Compute the solution change (absolute and relative changes).
    abs_schg   = norm(x_next - x_cur, 2);
    rel_schg   = abs_schg/max(1, norm(x_cur, 2));
    
    % Compute the feasibility gap (absotule and relative).
    abs_pfeas  = gy_val;
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
            x_cur         = x_next;
            break;
        end
    end
    % Update the second smoothness parameter beta.
    beta       = (1 - tau_next)*beta;

    % Perform Restarting if required.
    if options.isRestart
        if mod(iter, options.nRestart) == 0
            x_hat     = x_next;
            Ax_hat    = Ax_next;
            y_dot     = y_hat_star;
            beta      = beta1;
            tau_next  = 1.0;
        end
    end
    
    % Assign variables to go to the next step.
    tau       = tau_next;
    x_cur     = x_next;
    Ax_cur    = Ax_next;
end
% End of the main loop.

% Finalization.
PAPA_finalization;

% Get the final solution.
optsol.x_opt  = x_cur;
optsol.y_opt  = y_hat_star;
optsol.fx_val = fx_val;

end

% ASGARD v.1.0 by Quoc Tran-Dinh (quoctd@email.unc.edu)
% Joint work with Olivier Fercoq, and Volkan Cevher.
% Copyright @ 2016 Department of Statistics and Operations Research (STOR)
%                The University of North Carolina at Chapel Hill (UNC)
% See the file LICENSE for full license information.