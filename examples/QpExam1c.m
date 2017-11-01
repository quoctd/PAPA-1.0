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
mu_Q    = 1; % mu_Q = 1 for strongly convex case, and mu_Q = 0, otherwise.

width_r = 1;
p2      = round(p/2) + 1;
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
rsFreq1                = 100;
rsFreq2                = 50;
rsFreq3                = 25;

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
options.Algorithm  = 'PAPA-RS';
time1              = tic;
options.nRestart   = rsFreq1;
[optsol1, output1] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt1             = optsol1.y_opt;
fx_val1            = 0.5*y_opt1'*Q*y_opt1 + q'*y_opt1;
rel_feas1          = norm( min(B*y_opt1 - lbx, 0), 2 ) + norm( max(B*y_opt1 - ubx, 0), 2 );
diff_y1            = norm(y_opt1 - y_cvx, 2)/nrm_y_cvx;
time1              = toc(time1);

%% Call the PAPA solver with restart.
options.Algorithm  = 'PAPA-RS';
options.nRestart   = rsFreq2;
time2              = tic;
[optsol2, output2] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt2             = optsol2.y_opt;
fx_val2            = 0.5*y_opt2'*Q*y_opt2 + q'*y_opt2;
rel_feas2          = norm( min(B*y_opt2 - lbx, 0), 2 ) + norm( max(B*y_opt2 - ubx, 0), 2 );
diff_y2            = norm(y_opt2 - y_cvx, 2)/nrm_y_cvx;
time2              = toc(time2);

%% Call the PAPA solver without restart.
options.Algorithm  = 'PAPA-RS';
options.nRestart   = rsFreq3;
time3              = tic;
[optsol3, output3] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt3             = optsol3.y_opt;
fx_val3            = 0.5*y_opt3'*Q*y_opt3 + q'*y_opt3;
rel_feas3          = norm( min(B*y_opt3 - lbx, 0), 2 ) + norm( max(B*y_opt3 - ubx, 0), 2 );
diff_y3            = norm(y_opt3 - y_cvx, 2)/nrm_y_cvx;
time3              = toc(time3);

%% Call the PAPA solver without restart.
options.Algorithm  = 'PALPA-RS';
options.nRestart   = rsFreq1;
time4              = tic;
[optsol4, output4] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt4             = optsol4.y_opt;
fx_val4            = 0.5*y_opt4'*Q*y_opt4 + q'*y_opt4;
rel_feas4          = norm( min(B*y_opt4 - lbx, 0), 2 ) + norm( max(B*y_opt4 - ubx, 0), 2 );
diff_y4            = norm(y_opt4 - y_cvx, 2)/nrm_y_cvx;
time4              = toc(time4);

%% Call the PAPA solver for strongly convex case and without restart.
options.Algorithm  = 'PALPA-RS';
options.nRestart   = rsFreq2;
time5              = tic;
[optsol5, output5] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt5             = optsol5.y_opt;
fx_val5            = 0.5*y_opt5'*Q*y_opt5 + q'*y_opt5;
rel_feas5          = norm( min(B*y_opt5 - lbx, 0), 2 ) + norm( max(B*y_opt5 - ubx, 0), 2 );
diff_y5            = norm(y_opt5 - y_cvx, 2)/nrm_y_cvx;
time5              = toc(time5);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'PALPA-RS';
options.nRestart   = rsFreq3;
time6              = tic;
[optsol6, output6] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt6             = optsol6.y_opt;
fx_val6            = 0.5*y_opt6'*Q*y_opt6 + q'*y_opt6;
rel_feas6          = norm( min(B*y_opt6 - lbx, 0), 2 ) + norm( max(B*y_opt6 - ubx, 0), 2 );
diff_y6            = norm(y_opt6 - y_cvx, 2)/nrm_y_cvx;
time6              = toc(time6);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'PAPA-SCVX-RS';
options.nRestart   = rsFreq1;
time7              = tic;
[optsol7, output7] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt7             = optsol7.y_opt;
fx_val7            = 0.5*y_opt7'*Q*y_opt7 + q'*y_opt7;
rel_feas7          = norm( min(B*y_opt7 - lbx, 0), 2 ) + norm( max(B*y_opt7 - ubx, 0), 2 );
diff_y7            = norm(y_opt7 - y_cvx, 2)/nrm_y_cvx;
time7              = toc(time7);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'PAPA-SCVX-RS';
options.nRestart   = rsFreq2;
time8              = tic;
[optsol8, output8] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt8             = optsol8.y_opt;
fx_val8            = 0.5*y_opt8'*Q*y_opt8 + q'*y_opt8;
rel_feas8          = norm( min(B*y_opt8 - lbx, 0), 2 ) + norm( max(B*y_opt8 - ubx, 0), 2 );
diff_y8            = norm(y_opt8 - y_cvx, 2)/nrm_y_cvx;
time8              = toc(time8);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'PAPA-SCVX-RS';
options.nRestart   = rsFreq3;
time9              = tic;
[optsol9, output9] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt9             = optsol9.y_opt;
fx_val9            = 0.5*y_opt9'*Q*y_opt9 + q'*y_opt9;
rel_feas9          = norm( min(B*y_opt9 - lbx, 0), 2 ) + norm( max(B*y_opt9 - ubx, 0), 2 );
diff_y9            = norm(y_opt9 - y_cvx, 2)/nrm_y_cvx;
time9              = toc(time9);

%% Redefine the user-defined function for ASGARD ...
fzFunc2            = fzFunc;
fzFunc2.usDefFunc  = @(x, y, varargin) ( ( norm( min(B*x(n+1:end) - lbx, 0), 2 ) + norm( max(B*x(n+1:end) - ubx, 0), 2 ) ) /max_bld );

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm    = 'ASGARD-RS';
options.nRestart     = rsFreq1;
time10               = tic;
[optsol10, output10] = constrOptAsgardSolver(fzFunc2, B, cb, options, y0);
y_opt10              = optsol10.y_opt;
fx_val10             = 0.5*y_opt10'*Q*y_opt10 + q'*y_opt10;
rel_feas10           = norm( min(B*y_opt10 - lbx, 0), 2 ) + norm( max(B*y_opt10 - ubx, 0), 2 );
diff_y10             = norm(y_opt10 - y_cvx, 2)/nrm_y_cvx;
time10               = toc(time10);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm    = 'ASGARD-RS';
time11               = tic;
options.nRestart     = rsFreq2;
[optsol11, output11] = constrOptAsgardSolver(fzFunc2, B, cb, options, y0);
y_opt11              = optsol11.y_opt;
fx_val11             = 0.5*y_opt11'*Q*y_opt11 + q'*y_opt11;
rel_feas11           = norm( min(B*y_opt11 - lbx, 0), 2 ) + norm( max(B*y_opt11 - ubx, 0), 2 );
diff_y11             = norm(y_opt11 - y_cvx, 2)/nrm_y_cvx;
time11               = toc(time11);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm    = 'ASGARD-RS';
options.nRestart     = rsFreq3;
time12               = tic;
[optsol12, output12] = constrOptAsgardSolver(fzFunc2, B, cb, options, y0);
y_opt12              = optsol12.y_opt;
fx_val12             = 0.5*y_opt12'*Q*y_opt12 + q'*y_opt12;
rel_feas12           = norm( min(B*y_opt12 - lbx, 0), 2 ) + norm( max(B*y_opt12 - ubx, 0), 2 );
diff_y12             = norm(y_opt12 - y_cvx, 2)/nrm_y_cvx;
time12               = toc(time12);

%% Plot the outputs.
fx_min   = fx_cvx;
myabs    = @(x)( abs(x) );
YMatrix1 = [myabs(output1.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output2.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output3.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output4.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output5.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output6.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output7.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output8.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output9.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output10.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output11.hist.fx_val - fx_min)/max(1, abs(fx_min)), ...
            myabs(output12.hist.fx_val - fx_min)/max(1, abs(fx_min)) ...
            ];
YMatrix2 = [output1.hist.usedef, output2.hist.usedef, output3.hist.usedef, ...
            output4.hist.usedef, output5.hist.usedef, output6.hist.usedef, ...
            output7.hist.usedef, output8.hist.usedef, output9.hist.usedef, ...
            output10.hist.usedef, output11.hist.usedef, output12.hist.usedef];

%% Generate a figure
QpCreateFigure1c( YMatrix1(:,[1,2,7,8,10,11]), YMatrix2(:,[1,2,7,8,10,11]) );

% PAPA v.1.0 by Quoc Tran-Dinh (quoctd@email.unc.edu)
% Copyright @ 2017 Department of Statistics and Operations Research (STOR)
%                The University of North Carolina at Chapel Hill (UNC)
% See the file LICENSE for full license information.
