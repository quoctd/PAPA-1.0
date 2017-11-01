function [optsol, output] = QpPapaSolver(Q, q, B, lbx, ubx, options, ...
                            y0, varargin)
% This implement PAPA to solve the following QP problem:
%        min_y 0.5*y'*Q*y + q'*y s.t. lbx <= B*y <= ubx.
% We can reformulate this problem as
%        min g(x) + h(y) s.t. -x + B*y = 0.
%  where h(y) = 0.5*y'*Q*y + q'*y and g(x) = delta_[lbx, ubx](x).
%
if nargin < 5,  error('At least 5 inputs are required!\n'); end
if nargin < 6,  options = PAPA_OptimSet([]);    end
if nargin < 7,  y0      = zeros(size(B, 2), 1); end
if isempty(y0), y0      = zeros(size(B, 2), 1); end

% Check if h is strongly convex.
if strcmpi( options.Algorithm, 'PAPA-SCVX' ) || strcmpi( options.Algorithm, 'PAPA-SCVX-RS' )
    objFunc.muhy       = min(eig(Q));
    if objFunc.muhy <= 0
        warning('The strong convexity parameter must be positive!\n'); 
        objFunc.muhy   = 0.5;
    end
end

% This allows us to reduce the compuational complexity.
[V, S] = eig(Q); dS = diag(S);
x0     = B*y0;

% Define problem.
px                    = size(B, 1);
py                    = length(q);
objFunc.nc            = px;
objFunc.px            = px;
objFunc.py            = py;
objFunc.gxProxOper    = @(x, gamma)( min(max(x, lbx), ubx) );
objFunc.hyProxOper    = @(y, gamma, varargin) ( V*((1./(1 + gamma*dS)).*(V'*(y - gamma*q))) );
objFunc.gxFunc        = @(x, gamma, varargin) ( 0 );
objFunc.hyFunc        = @(y, gamma, varargin) ( 0.5*y'*Q*y + q'*y );
linConstr.Aoper       = @(x, varargin) ( -x );
linConstr.AToper      = @(u, varargin) ( -u );
linConstr.Boper       = @(y, varargin) ( B*y );
linConstr.BToper      = @(v, varargin) ( B'*v );
linConstr.dlStarProx  = @(u, varargin) ( u );
linConstr.dlStarFunc  = @(u, varargin) ( norm(u, 2) );

% Set other parameters.
LA_bar                = 1;
LB_bar                = PAPA_l2NormEval(objFunc.py, linConstr.Boper, ...
                        linConstr.BToper, options.PwMaxIters, options.PwRelTol);
linConstr.LA_bar      = LA_bar;
linConstr.LB_bar      = LB_bar;
    
%% Call the solver ...
if strcmpi(options.Algorithm, 'PAPA')
    
    % The PAPA solver without restart
    options_ol           = options;
    options_ol.isRestart = 0;
    options_ol.beta1     = sqrt(LB_bar);
    [optsol, output]     = nscvxPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PAPA-RS')
    
    % The PAPA solver with restart
    options_ol           = options;
    options_ol.isRestart = 1;
    options_ol.beta1     = sqrt(LB_bar);
    [optsol, output]     = nscvxPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PALPA') 
    
   % The linearized PAPA (PALPA) solver without restart
    options_ol           = options;
    options_ol.isRestart = 0;
    options_ol.beta1     = sqrt(LA_bar + LB_bar);
    [optsol, output]     = nscvxLinearizedPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PALPA-RS') 
    
   % The linearized PAPA (PALPA) solver with restart
    options_ol           = options;
    options_ol.isRestart = 1;
    options_ol.beta1     = sqrt(LA_bar + LB_bar);
    [optsol, output]     = nscvxLinearizedPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PAPA-SCVX') 
    
   % The PAPA solver for strongly convex case and without restart
    options_ol           = options;
    options_ol.isRestart = 0;
    options_ol.beta1     = 2.0*LB_bar/objFunc.muhy;
    [optsol, output]     = scvxPapaSolver(objFunc, linConstr, x0, y0, options_ol);
    
elseif strcmpi(options.Algorithm, 'PAPA-SCVX-RS') 
   
   % The PAPA solver for strongly convex case with restart
    options_ol           = options;
    options_ol.isRestart = 1;
    options_ol.beta1     = 2.0*LB_bar/objFunc.muhy;
    [optsol, output]     = scvxPapaSolver(objFunc, linConstr, x0, y0, options_ol);

else
    warning('Unknown solver!\n');
    optsol = []; output = [];
end

% PAPA v.1.0 by Quoc Tran-Dinh (quoctd@email.unc.edu)
% Copyright @ 2017 Department of Statistics and Operations Research (STOR)
%                The University of North Carolina at Chapel Hill (UNC)
% See the file LICENSE for full license information.