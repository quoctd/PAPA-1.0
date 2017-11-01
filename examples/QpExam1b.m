% Generate the input data.
%% Example: Quadratic programming:
%
%          min 0.5*y'*Q*y + q'*y s.t. lbx <= B*y <= ubx.
%
%  Here, Q in R^{pxp} is symmetric positive semidefinite, B in R^{nxp},
%  q in R^p, and lbx, ubx in R^n.
%

p       = 200;
n       = 200;
mu_Q    = 1;  % mu_Q = 1 for strongly convex case, and mu_Q = 0, otherwise.

width_r = 1;
p2      = round(0.5*p) + 1;
R       = randn(p, p2)/sqrt(p2);
Q_nscvx = R*R';
q       = randn(p,1);

B       = randn(n, p);
B       = B/sqrt(n);
x_org   = randn(p,1);
lbx     = B*x_org - width_r*rand(n,1);
ubx     = B*x_org + width_r*rand(n,1);
cb      = zeros(n, 1);
max_bld = max([norm(lbx,2), norm(ubx,2), 1]);

% Generate measurement.
Q = Q_nscvx + mu_Q*eye(p);

%% Define the objective function and its proximal operators.
[V, S] = eig(Q); 
dS     = diag(S);
fzFunc.gxProxOper   = @(x, gamma)( min(max(x, lbx), ubx) );
fzFunc.gxFunc       = @(x, gamma, varargin) ( 0 );
fzFunc.hyProxOper   = @(y, gamma, varargin) ( V*((1./(1 + gamma*dS)).*(V'*(y - gamma*q))) );
fzFunc.hyFunc       = @(y, gamma, varargin) ( 0.5*y'*Q*y + q'*y );
fzFunc.usDefFunc    = @(x, y, varargin) ( ( norm( min(B*y - lbx, 0), 2 ) + norm( max(B*y - ubx, 0), 2 ) ) /max_bld );

%% Define the convexity parameter.
muhy   = min(dS);
if muhy <= 0, muhy = 0.1; end
fzFunc.muhy = muhy;

%% Generate a starting point.
y0 = zeros(p, 1);

% Set the optional parameters.
options                = PAPA_OptimSet([]);
options.isStoppingCond = 0;
options.saveHistMode   = 4;
options.MaxIters       = 1000;
rsFreq                 = 50;
rsFreq2                = 100;

%% Solve the problem by CVX.
if exist('cvx_begin.m')
    time_cvx = tic;
    cvx_solver mosek;
    cvx_begin;
        variable y_cvx(p);
        minimize( 0.5*quad_form(y_cvx, Q) + q'*y_cvx );
        subject to
            lbx <= B*y_cvx; B*y_cvx <= ubx;
    cvx_end;
    time_cvx     = toc(time_cvx);
    fx_cvx       = 0.5*y_cvx'*Q*y_cvx + q'*y_cvx;
    rel_feas_cvx = norm( min(B*y_cvx    - lbx, 0), 2 ) + norm( max(B*y_cvx    - ubx, 0), 2 );
    nrm_y_cvx    = max(1, norm(y_cvx, 2));
else
    % Call Quadprog.
    time_cvx        = tic;
    opts            = optimset('Algorithm', 'interior-point-convex', 'TolCon', 1e-14, 'TolConSQP', 1e-12, 'Display', 'iter');
    [y_cvx, fx_cvx] = quadprog(Q, q, [B; -B], [ubx; -lbx], [], [], [], [], [], opts);
    time_cvx        = toc(time_cvx);
    fx_cvx          = 0.5*y_cvx'*Q*y_cvx + q'*y_cvx;
    rel_feas_cvx    = norm( min(B*y_cvx    - lbx, 0), 2 ) + norm( max(B*y_cvx    - ubx, 0), 2 );
    nrm_y_cvx       = max(1, norm(y_cvx, 2));
end

%% Call the PAPA solver without restart.
options.Algorithm  = 'PAPA';
time1              = tic;
[optsol1, output1] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt1             = optsol1.y_opt;
fx_val1            = 0.5*y_opt1'*Q*y_opt1 + q'*y_opt1;
rel_feas1          = norm( min(B*y_opt1 - lbx, 0), 2 ) + norm( max(B*y_opt1 - ubx, 0), 2 );
diff_y1            = norm(y_opt1 - y_cvx, 2)/nrm_y_cvx;
time1              = toc(time1);

%% Call the PAPA solver with restart.
options.Algorithm  = 'PAPA-RS';
options.nRestart   = rsFreq;
time2              = tic;
[optsol2, output2] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt2             = optsol2.y_opt;
fx_val2            = 0.5*y_opt2'*Q*y_opt2 + q'*y_opt2;
rel_feas2          = norm( min(B*y_opt2 - lbx, 0), 2 ) + norm( max(B*y_opt2 - ubx, 0), 2 );
diff_y2            = norm(y_opt2 - y_cvx, 2)/nrm_y_cvx;
time2              = toc(time2);

%% Call the PAPA solver without restart.
options.Algorithm  = 'PALPA';
time3              = tic;
[optsol3, output3] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt3             = optsol3.y_opt;
fx_val3            = 0.5*y_opt3'*Q*y_opt3 + q'*y_opt3;
rel_feas3          = norm( min(B*y_opt3 - lbx, 0), 2 ) + norm( max(B*y_opt3 - ubx, 0), 2 );
diff_y3            = norm(y_opt3 - y_cvx, 2)/nrm_y_cvx;
time3              = toc(time3);

%% Call the PAPA solver without restart.
options.Algorithm  = 'PALPA-RS';
options.nRestart   = rsFreq;
time4              = tic;
[optsol4, output4] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt4             = optsol4.y_opt;
fx_val4            = 0.5*y_opt4'*Q*y_opt4 + q'*y_opt4;
rel_feas4          = norm( min(B*y_opt4 - lbx, 0), 2 ) + norm( max(B*y_opt4 - ubx, 0), 2 );
diff_y4            = norm(y_opt4 - y_cvx, 2)/nrm_y_cvx;
time4              = toc(time4);

%% Call the PAPA solver for strongly convex case and without restart.
options.Algorithm  = 'PAPA-SCVX';
time5              = tic;
[optsol5, output5] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt5             = optsol5.y_opt;
fx_val5            = 0.5*y_opt5'*Q*y_opt5 + q'*y_opt5;
rel_feas5          = norm( min(B*y_opt5 - lbx, 0), 2 ) + norm( max(B*y_opt5 - ubx, 0), 2 );
diff_y5            = norm(y_opt5 - y_cvx, 2)/nrm_y_cvx;
time5              = toc(time5);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'PAPA-SCVX-RS';
options.nRestart   = rsFreq2;
time6              = tic;
[optsol6, output6] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt6             = optsol6.y_opt;
fx_val6            = 0.5*y_opt6'*Q*y_opt6 + q'*y_opt6;
rel_feas6          = norm( min(B*y_opt6 - lbx, 0), 2 ) + norm( max(B*y_opt6 - ubx, 0), 2 );
diff_y6            = norm(y_opt6 - y_cvx, 2)/nrm_y_cvx;
time6              = toc(time6);

%% Redefine the user-defined function for ASGARD ...
fzFunc2            = fzFunc;
fzFunc2.usDefFunc  = @(x, y, varargin) ( ( norm( min(B*x(n+1:end) - lbx, 0), 2 ) + norm( max(B*x(n+1:end) - ubx, 0), 2 ) ) /max_bld );

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'ASGARD';
time7              = tic;
[optsol7, output7] = constrOptAsgardSolver(fzFunc2, B, cb, options, y0);
y_opt7             = optsol7.y_opt;
fx_val7            = 0.5*y_opt7'*Q*y_opt7 + q'*y_opt7;
rel_feas7          = norm( min(B*y_opt7 - lbx, 0), 2 ) + norm( max(B*y_opt7 - ubx, 0), 2 );
diff_y7            = norm(y_opt7 - y_cvx, 2)/nrm_y_cvx;
time7              = toc(time7);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'ASGARD-RS';
options.nRestart   = rsFreq;
time8              = tic;
[optsol8, output8] = constrOptAsgardSolver(fzFunc2, B, cb, options, y0);
y_opt8             = optsol8.y_opt;
fx_val8            = 0.5*y_opt8'*Q*y_opt8 + q'*y_opt8;
rel_feas8          = norm( min(B*y_opt8 - lbx, 0), 2 ) + norm( max(B*y_opt8 - ubx, 0), 2 );
diff_y8            = norm(y_opt8 - y_cvx, 2)/nrm_y_cvx;
time8              = toc(time8);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'CP';
time9              = tic;
fzFunc3            = fzFunc;
fzFunc3.usDefFunc  = @(x, y, varargin) ( ( norm( min(B*x - lbx, 0), 2 ) + norm( max(B*x - ubx, 0), 2 ) ) /max_bld );
[optsol9, output9] = constrOptAsgardSolver(fzFunc3, B, cb, options, y0);
y_opt9             = optsol9.x_opt;
fx_val9            = 0.5*y_opt9'*Q*y_opt9 + q'*y_opt9;
rel_feas9          = norm( min(B*y_opt9 - lbx, 0), 2 ) + norm( max(B*y_opt9 - ubx, 0), 2 );
diff_y9            = norm(y_opt9 - y_cvx, 2)/nrm_y_cvx;
time9              = toc(time9);

%% Print the final results.
fprintf(['The objective values:       \n', ...
         '  +PAPA         = %3.20f\n  +PAPA-RS      = %3.20f\n', ...
         '  +PALPA        = %3.20f\n  +PALPA-RS     = %3.20f\n  +PAPA-SCVX    = %3.20f\n', ...
         '  +PAPA-SCVX-RS = %3.20f\n  +ASGARD       = %3.20f\n  +ASGARD-RS    = %3.20f\n', ...
         '  +CP           = %3.20f\n  +CVX          = %3.20f\n'], ...
         fx_val1, fx_val2, fx_val2, fx_val4, fx_val5, fx_val6, ...
         fx_val7, fx_val8, fx_val9, fx_cvx);
fprintf(['The absolute feasibility gap:       \n', ...
         '  +PAPA         = %3.20f\n  +PAPA-RS      = %3.20f\n', ...
         '  +PALPA        = %3.20f\n  +PALPA-RS     = %3.20f\n  +PAPA-SCVX    = %3.20f\n', ...
         '  +PAPA-SCVX-RS = %3.20f\n  +ASGARD       = %3.20f\n  +ASGARD-RS    = %3.20f\n', ...
         '  +CP           = %3.20f\n  +CVX          = %3.20f\n'], ...
         rel_feas1, rel_feas2, rel_feas3, rel_feas4, rel_feas5, rel_feas6, ...
         rel_feas7, rel_feas8, rel_feas9, rel_feas_cvx);
fprintf(['The solution differences:\n', ...
         '  +PAPA         = %3.20f\n  +PAPA-RS      = %3.20f\n', ...
         '  +PALPA        = %3.20f\n  +PALPA-RS     = %3.20f\n  +PAPA-SCVX    = %3.20f\n', ...
         '  +PAPA-SCVX-RS = %3.20f\n  +ASGARD       = %3.20f\n  +ASGARD-RS    = %3.20f\n', ...
         '  +CP           = %3.20f\n'], ...   
         diff_y1, diff_y2, diff_y3, diff_y4, diff_y5, diff_y6, diff_y7, diff_y8, diff_y9);
fprintf(['The solution time in second:\n', ...
         '  +PAPA         = %3.20f\n  +PAPA-RS      = %3.20f\n', ...
         '  +PALPA        = %3.20f\n  +PALPA-RS     = %3.20f\n  +PAPA-SCVX    = %3.20f\n', ...
         '  +PAPA-SCVX-RS = %3.20f\n  +ASGARD       = %3.20f\n  +ASGARD-RS    = %3.20f\n', ...
         '  +CP           = %3.20f\n  +CVX          = %3.20f\n'], ...
         time1, time2, time3, time4, time5, time6, time7, time8, time9, time_cvx);

%% Plot the outputs.
fx_min   = fx_cvx;
myabs    = @(x)( abs(x) );
YMatrix1 = [myabs(output1.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output2.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output5.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output6.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output7.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output8.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output9.hist.fx_val - fx_min)/max(1, abs(fx_min))];
YMatrix2 = [output1.hist.usedef, output2.hist.usedef, output5.hist.usedef, output6.hist.usedef, ...
            output7.hist.usedef, output8.hist.usedef, output9.hist.usedef];


QpCreateFigure1b( YMatrix1, YMatrix2 );

% PAPA v.1.0 by Quoc Tran-Dinh (quoctd@email.unc.edu)
% Copyright @ 2017 Department of Statistics and Operations Research (STOR)
%                The University of North Carolina at Chapel Hill (UNC)
% See the file LICENSE for full license information.


