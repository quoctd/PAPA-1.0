function [optsol, output] = constrOptAsgardSolver(fzFunc, B, c, options, ...
                                                  y0, varargin)
% This implement Asgard to solve the following constrained convex problem:
%
%        min_{x, y} F(z) = g(x) + h(y) s.t. -x + B*y = c.
% Here:
%    f and g are convex and have tractable proximity operator.
%
if nargin < 3, error('At least 6 inputs are required!\n'); end
if nargin < 4, options = PAPA_OptimSet([]); end
if nargin < 5,    y0   = zeros(size(B, 2), 1); end
if isempty(y0),   y0   = zeros(size(B, 2), 1); end

% Generate an initial point for x0.
x0     = [B*y0 - c; y0];

% Define problem.
px                    = size(B, 1);
py                    = size(B, 2);
objFunc.nx            = px + py;
objFunc.fxProxOper    = @(x, gamma, varargin) ( [fzFunc.gxProxOper(x(1:px), gamma); fzFunc.hyProxOper(x(px+1:end), gamma)] );
objFunc.fxFunc        = @(x, gamma, varargin) ( fzFunc.gxFunc( x(1:px) ) + fzFunc.hyFunc( x(px+1:end) ) );
if isfield(fzFunc, 'usDefFunc'), objFunc.usDefFunc = fzFunc.usDefFunc; end
if isfield(fzFunc, 'muhy')
    if fzFunc.muhy <= 0 
        warning('The strong convexity parameter must be positive'); 
    end
    objFunc.muhy      = fzFunc.muhy;
end    
linConstr.Aoper       = @(x, varargin) ( -x(1:px) + B*x(px+1:end) );
linConstr.AToper      = @(y, varargin) ( [-y; B'*y]);
linConstr.cb          = c;
linConstr.gyStarProx  = @(u, gamma)( u - gamma*c  );
linConstr.gyStarFunc  = @(u, varargin) ( norm(u - c, 2) );

% Set other parameters.
LA_bar                  = PAPA_l2NormEval(objFunc.nx, linConstr.Aoper, ...
                          linConstr.AToper, options.PwMaxIters, options.PwRelTol);
options.LA_bar          = LA_bar;
options.isStoppingCond  = false;
trade_off1              = 1.0;
    
%% Call the solver ...
if strcmpi(options.Algorithm, 'ASGARD')
    
    % One Loop Asgard without restart
    options_ol             = options;
    options_ol.isRestart   = 0;
    options_ol.beta1       = trade_off1*sqrt(LA_bar);

    [optsolxy, output]     = copAsgardSolver(objFunc, linConstr, x0, options_ol);
    optsol.x_opt           = optsolxy.x_opt(1:px);
    optsol.y_opt           = optsolxy.x_opt(px+1:end);
    
elseif strcmpi(options.Algorithm, 'ASGARD-RS') 
    
    % One Loop Asgard restart
    options_olr           = options;
    options_olr.isRestart = 1;
    options_olr.nRestart  = 25;
    options_olr.beta1     = trade_off1*sqrt(LA_bar);
    [optsolxy, output]    = copAsgardSolver(objFunc, linConstr, x0, options_olr);
    optsol.x_opt          = optsolxy.x_opt(1:px);
    optsol.y_opt          = optsolxy.x_opt(px+1:end);
   
elseif strcmpi(options.Algorithm, 'CP') 
    
    % Define problem.
    objFunc.nx            = py;
    objFunc.fxProxOper    = @(x, gamma, varargin) ( fzFunc.gxProxOper(x - c, gamma, varargin{:}) + c );
    objFunc.gxProxOper    = @(u, gamma, varargin) ( fzFunc.hyProxOper(u, gamma, varargin{:}) );
    objFunc.fxFunc        = @(x, varargin) ( fzFunc.gxFunc( x, varargin{:}) );
    objFunc.gxFunc        = @(u, varargin) ( fzFunc.hyFunc(u, varargin{:}) );
    
    linConstr.Aoper       = @(x, varargin) ( B*x  );
    linConstr.AToper      = @(x, varargin) ( B'*x );
    if isfield(fzFunc, 'usDefFunc'), objFunc.usDefFunc = fzFunc.usDefFunc; end
    if isfield(fzFunc, 'muhy')
        if fzFunc.muhy <= 0 
            warning('The strong convexity parameter must be positive'); 
        end
        objFunc.muhy      = fzFunc.muhy;
    end
    
    % Chambolle-Pock method.
    if fzFunc.muhy > 0
        [optsolxy, output] = scvxChambollePockSolver(objFunc, linConstr, y0, options);
    else
        [optsolxy, output] = nscvxChambollePockSolver(objFunc, linConstr, y0, options);
    end
    optsol.x_opt       = optsolxy.x_opt;
    optsol.y_opt       = optsolxy.y_opt;
    
else
    
    % No solver is selected.
    warning('No solver is selected!\n');
    optsol = []; output = [];

end

% End of the implementation.