function [optsol, output] = constrOptPapaSolver(fzFunc, B, c, options, ...
                                                y0, varargin)
% This code implements PAPA to solve the following constrained convex problem:
%
%        min_{x, y} F(z) = g(x) + h(y) s.t. -x + B*y = c.
%
% Here:
%    f and g are convex and have tractable proximity operator.
%
%
if nargin < 3, error('At least 6 inputs are required!\n'); end
if nargin < 4, options = PAPA_OptimSet([]); end
if nargin < 5,    y0   = zeros(size(B, 2), 1); end
if isempty(y0),   y0   = zeros(size(B, 2), 1); end

% Check if h is strongly convex.
if strcmpi( options.Algorithm, 'PAPA-SCVX' ) || strcmpi( options.Algorithm, 'PAPA-SCVX-RS' )
    if ~isfield(fzFunc, 'muhy')
        warning('A strong convexity parameter is missing!'); 
    end
    if fzFunc.muhy <= 0 
        warning('The strong convexity parameter must be positive'); 
    end
    objFunc.muhy      = fzFunc.muhy;
end

% Generate an initial point for x0.
x0     = B*y0 - c;

% Define the problem objective function.
px                    = size(B, 1);
py                    = size(B, 2);
objFunc.nc            = px;
objFunc.px            = px;
objFunc.py            = py;
objFunc.gxProxOper    = @(x, gamma, varargin) ( fzFunc.gxProxOper(x - c, gamma, varargin{:}) );
objFunc.hyProxOper    = @(y, gamma, varargin) ( fzFunc.hyProxOper(y, gamma, varargin{:}) );
objFunc.gxFunc        = @(x, varargin) ( fzFunc.gxFunc(x, varargin{:}) );
objFunc.hyFunc        = @(y, varargin) ( fzFunc.hyFunc(y, varargin{:}) );
if isfield(fzFunc, 'usDefFunc'), objFunc.usDefFunc = fzFunc.usDefFunc; end

% Define linear constraints.
linConstr.Aoper       = @(x, varargin) ( -x );
linConstr.AToper      = @(u, varargin) ( -u );
linConstr.Boper       = @(y, varargin) ( B*y );
linConstr.BToper      = @(v, varargin) ( B'*v );
linConstr.dlStarProx  = @(u, gamma, varargin) ( u - gamma*c );
linConstr.dlStarFunc  = @(u, varargin) ( norm(u - c, 2) );

% Set other parameters.
LA_bar                  = 1;
LB_bar                  = PAPA_l2NormEval(objFunc.py, linConstr.Boper, ...
                          linConstr.BToper, options.PwMaxIters, options.PwRelTol);
linConstr.LA_bar        = LA_bar;
linConstr.LB_bar        = LB_bar;
trade_off1              = 1;
trade_off2              = 1;
    
%% Call the solver ...
if strcmpi(options.Algorithm, 'PAPA')
    
    % The PAPA solver without restart
    options_ol           = options;
    options_ol.isRestart = 0;
    options_ol.beta1     = trade_off1*sqrt(LB_bar);
    
    % Call the solver.
    [optsol, output]     = nscvxPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PAPA-RS')
    
    % The PAPA solver with restart
    options_ol           = options;
    options_ol.isRestart = 1;
    options_ol.beta1     = trade_off1*sqrt(LB_bar);
    
    % Call the solver.
    [optsol, output]     = nscvxPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PALPA') 
    
    % The linearized PAPA (PALPA) solver without restart
    objFunc.gxProxOper    = @(x, gamma, varargin) ( fzFunc.gxProxOper(x, gamma, varargin{:}) );
   
    % Set optional parameters and call the solver.
    options_ol           = options;
    options_ol.isRestart = 0;
    options_ol.beta1     = trade_off2*sqrt(LA_bar + LB_bar);
    
    % Call the solver.
    [optsol, output]     = nscvxLinearizedPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PALPA-RS') 
    
    % The linearized PAPA (PALPA) solver with restart
    objFunc.gxProxOper    = @(x, gamma, varargin) ( fzFunc.gxProxOper(x, gamma, varargin{:}) );
    
    % Set optional parameters and call the solver.
    options_ol           = options;
    options_ol.isRestart = 1;
    options_ol.beta1     = trade_off2*sqrt(LA_bar + LB_bar);
    
    % Call the solver.
    [optsol, output]     = nscvxLinearizedPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PAPA-SCVX') 
    
   % The PAPA solver for strongly convex case and without restart
    options_ol           = options;
    options_ol.isRestart = 0;
    options_ol.beta1     = 2.0*LB_bar/max(objFunc.muhy, options.lbScvxParam);
    
    % Call the solver.
    [optsol, output]     = scvxPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PAPA-SCVX-RS') 
   
   % The PAPA solver for strongly convex case with restart
    options_ol           = options;
    options_ol.isRestart = 1;
    options_ol.beta1     = 2.0*LB_bar/max(objFunc.muhy, options.lbScvxParam);
    
    % Call the solver.
    [optsol, output]     = scvxPapaSolver(objFunc, linConstr, x0, y0, options_ol);

else
    warrning('Unknown method!\n');
    optsol = []; output = [];
end

% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.