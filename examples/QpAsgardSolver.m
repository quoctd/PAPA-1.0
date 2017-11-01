function [optsol, output] = QpAsgardSolver(Q, q, B, c, lbx, ubx, ...
                                           options, y0, varargin)
% This implement Asgard to solve the following QP problem:
%        min_y 0.5*y'*Q*y + q'*y s.t. lbx <= B*y <= ubx.
% We can reformulate this problem as
%        min f(y) + g(B*y),
%  where f(y) = 0.5*y'*Q*y + q'*y and g(u) = delta_[lbx, ubx](u).
%
if nargin < 6, error('At least 6 inputs are required!\n'); end
if nargin < 7, options = PAPA_OptimSet([]); end
if nargin < 8,    y0   = zeros(size(B, 2), 1); end
if isempty(y0),   y0   = zeros(size(B, 2), 1); end

% This allows us to reduce the compuational complexity.
[V, S] = eig(Q); dS = diag(S);

% Define problem.
nx                    = length(q);
objFunc.nx            = nx;
objFunc.fxProxOper    = @(x, gamma, varargin) ( V*((1./(1 + gamma*dS)).*(V'*(x - gamma*q))) );
objFunc.fxFunc        = @(x, gamma, varargin) ( 0.5*x'*Q*x + q'*x );
linConstr.Aoper       = @(x, varargin) B*x;
linConstr.AToper      = @(x, varargin) B'*x;
linConstr.cb          = c;
linConstr.gyStarProx  = @(x, gamma)( x - gamma * min(max(x/gamma, lbx + c), ubx + c) );
linConstr.gyStarFunc  = @(u, varargin) ( norm(max(u - c - ubx, 0), 2) + norm(min(u - c - lbx, 0), 2) );


% Set other parameters.
LA_bar                  = ASGARD_l2NormEval(objFunc.nx, linConstr.Aoper, ...
                          linConstr.AToper, options.PwMaxIters, options.PwRelTol);
options.LA_bar          = LA_bar;
options.isStoppingCond  = false;
    
%% Call the solver ...
if strcmpi(options.Algorithm, '2LASGARD')
    
    % Two Loop Asgard
    alpha            = 1.1;
    m_0              = max(10, round(alpha/(alpha-1)) + 1);
    outer_dl         = 100;
    options.num_eps  = outer_dl;
    options_dl       = options;
    options_dl.beta1 = 0.05*sqrt(LA_bar);
    options_dl.m_0   = m_0;
    options_dl.alpha = alpha;

    [optsol, output] = TwoLoopAsgardSolver(objFunc, linConstr, y0, options_dl);
    
elseif strcmpi(options.Algorithm, '1LASGARD')
    
    % One Loop Asgard
    options_ol           = options;
    options_ol.isRestart = 0;
    options_ol.beta1     = 1*sqrt(LA_bar);

    [optsol, output]     = copAsgardSolver(objFunc, linConstr, y0, options_ol);
    
elseif strcmpi(options.Algorithm, '1LRASGARD') 
    
    % One Loop Asgard restart
    options_olr           = options;
    options_olr.isRestart = 1;
    options_olr.nRestart  = 25;
    options_olr.beta1     = 1*sqrt(LA_bar);
    [optsol, output]      = copAsgardSolver(objFunc, linConstr, y0, options_olr);
    
else
    
    % Chambolle-Pock method.
    [optsol, output] = ChambollePockSolver(objFunc, linConstr, y0, options);
    
end

% End of the implementation.