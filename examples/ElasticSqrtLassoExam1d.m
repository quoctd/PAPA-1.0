% Example: Solve the Sqrt-LASSO when setting reg1 = 0.
% To run this code, CVX needs to be installed. Otherwise, it takes zero as
% the default objective value.
% See: http://cvxr.com/

% Generate the input data.
scale           = 10; % Change this to get different problem size
p               = scale*100;
n               = scale*35; 
s               = scale*10; 

% Generate matrix B.
cor_tau = 0.0;
if cor_tau > 0
  var0 = (1 - cor_tau)^2 / (1 - cor_tau^2); %initial variance
  B = zeros(n, p);
  B(:,1) = sqrt(var0)*randn(n, 1);
  for kk = 2:p
    B(:,kk) = cor_tau*B(:,kk-1) + (1 - cor_tau)*(randn(n,1));
  end
else
    B   = randn(n, p);
end
B  = B/sqrt(n);

% Generate vector c.
x_org           = zeros(p, 1);
T               = randsample(n, s);
x_org(T)        = randn(s, 1);

% Generate problem noise and measurement
noise           = randn(n, 1);
reg1            = 0.0;
reg2            = 1e-1;
muhy            = 0.0;

%% TWO RUN.
for run = 1:2

% Generate measurement.
if run == 1, sigma = 0; else sigma = 1e-3; end
cb = B*x_org + sigma*noise;

%% Solving the problem
%         min |B*y - c|_2 + 0.5*r1*|y|^2 + r2*|y|_1.
% Reformulate this as
%         min |x|_2 + 0.5*r1*|y|^2 + r2*|y|_1 s.t. -x + B*y = c.
%
fzFunc.gxProxOper   = @(x, gamma, varargin) ( max( 1 - gamma/norm(x, 2), 0)*x );
fzFunc.gxFunc       = @(x, varargin) ( norm(x, 2) );
soft_thres_hold     = @(y, gamma) sign(y).*max(abs(y) - gamma, 0);
fzFunc.hyProxOper   = @(y, gamma, varargin) ( soft_thres_hold( y/(1 + reg1*gamma), reg2*gamma/(1 + reg1*gamma)) );
fzFunc.hyFunc       = @(y, varargin) ( 0.5*reg1*norm(y, 2)^2 + reg2*norm(y, 1) );
fzFunc.muhy         = muhy;
fzFunc.usDefFunc    = @(x, y, varargin) norm(B*y - cb, 2) + 0.5*reg1*norm(y, 2).^2 + reg2*norm(y, 1);

%% Generate a starting point.
y0                  = zeros(p, 1);

% Set the optional parameters.
options                = PAPA_OptimSet([]);
options.isStoppingCond = 0;
options.saveHistMode   = 4;
options.MaxIters       = 1000;
rsFreq                 = 50;
rsFreq2                = 100;

%% Call the CVX-Mosek solver.
if exist('cvx_begin.m')
    time_cvx = tic;
    cvx_solver mosek; % If mosek is not installed, comment this line.
    cvx_precision best;
    cvx_begin
        variable y_cvx(p);
        variable x_cvx(n);
        minimize( norm(x_cvx, 2) + 0.5*reg1*(y_cvx'*y_cvx) + reg2 * norm(y_cvx, 1));
        subject to
            x_cvx == B*y_cvx - cb;
    cvx_end
    fx_cvx    = norm(B*y_cvx - cb, 2) + 0.5*reg1*(y_cvx'*y_cvx) + reg2 * norm(y_cvx, 1);
    time_cvx  = toc(time_cvx);
    nrm_y_cvx = max(1, norm(y_cvx, 2)); 
else
    fx_cvx = 0; time_cvx = 0; nrm_y_cvx = nan;
end

%% Call the PAPA solver for strongly convex case and without restart.
if reg1 == 0, options.lbScvxParam = 1.0; end
options.Algorithm  = 'PAPA-SCVX';
time1              = tic;
[optsol1, output1] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt1             = optsol1.y_opt;
fx_val1            = norm(B*y_opt1 - cb, 2) + 0.5*reg1*norm(y_opt1, 2)^2 + reg2*norm(y_opt1, 1);
diff_y1            = norm(y_opt1 - y_cvx, 2)/nrm_y_cvx;
time1              = toc(time1);

%% Call the PAPA solver for strongly convex case and without restart.
if reg1 == 0, options.lbScvxParam = 1.0e-1; end
options.Algorithm  = 'PAPA-SCVX';
time2              = tic;
[optsol2, output2] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt2             = optsol2.y_opt;
fx_val2            = norm(B*y_opt2 - cb, 2) + 0.5*reg1*norm(y_opt2, 2)^2 + reg2*norm(y_opt2, 1);
diff_y2            = norm(y_opt2 - y_cvx, 2)/nrm_y_cvx;
time2              = toc(time2);

%% Call the PAPA solver for strongly convex case and without restart.
if reg1 == 0, options.lbScvxParam = 1.0e-2; end
options.Algorithm  = 'PAPA-SCVX';
time3              = tic;
[optsol3, output3] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt3             = optsol3.y_opt;
fx_val3            = norm(B*y_opt3 - cb, 2) + 0.5*reg1*norm(y_opt3, 2)^2 + reg2*norm(y_opt3, 1);
diff_y3            = norm(y_opt3 - y_cvx, 2)/nrm_y_cvx;
time3              = toc(time3);

%% Call the PAPA solver for strongly convex case and with restart.
if reg1 == 0, options.lbScvxParam = 1.0; end
options.Algorithm  = 'PAPA-SCVX-RS';
options.nRestart   = rsFreq2;
time4              = tic;
[optsol4, output4] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt4             = optsol4.y_opt;
fx_val4            = norm(B*y_opt4 - cb, 2) + 0.5*reg1*norm(y_opt4, 2)^2 + reg2*norm(y_opt4, 1);
diff_y4            = norm(y_opt4 - y_cvx, 2)/nrm_y_cvx;
time4              = toc(time4);

%% Call the PAPA solver for strongly convex case and with restart.
if reg1 == 0, options.lbScvxParam = 1.0e-1; end
options.Algorithm  = 'PAPA-SCVX-RS';
options.nRestart   = rsFreq2;
time5              = tic;
[optsol5, output5] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt5             = optsol5.y_opt;
fx_val5            = norm(B*y_opt5 - cb, 2) + 0.5*reg1*norm(y_opt5, 2)^2 + reg2*norm(y_opt5, 1);
diff_y5            = norm(y_opt5 - y_cvx, 2)/nrm_y_cvx;
time5              = toc(time5);

%% Call the PAPA solver for strongly convex case and with restart.
if reg1 == 0, options.lbScvxParam = 1.0e-2; end
options.Algorithm  = 'PAPA-SCVX-RS';
options.nRestart   = rsFreq2;
time6              = tic;
[optsol6, output6] = constrOptPapaSolver(fzFunc, B, cb, options, y0);
y_opt6             = optsol6.y_opt;
fx_val6            = norm(B*y_opt6 - cb, 2) + 0.5*reg1*norm(y_opt6, 2)^2 + reg2*norm(y_opt6, 1);
diff_y6            = norm(y_opt6 - y_cvx, 2)/nrm_y_cvx;
time6              = toc(time6);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'CP';
time7              = tic;
fzFunc.muhy        = 1.0;
fzFunc.usDefFunc    = @(x, y, varargin) norm(B*x - cb, 2) + 0.5*reg1*norm(x, 2).^2 + reg2*norm(x, 1);
[optsol7, output7] = constrOptAsgardSolver(fzFunc, B, cb, options, y0);
y_opt7             = optsol7.x_opt;
fx_val7            = norm(B*y_opt7 - cb, 2) + 0.5*reg1*norm(y_opt7, 2)^2 + reg2*norm(y_opt7, 1);
diff_y7            = norm(y_opt7 - y_cvx, 2)/nrm_y_cvx;
time7              = toc(time7);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'CP';
time8              = tic;
fzFunc.usDefFunc   = @(x, y, varargin) norm(B*x - cb, 2) + 0.5*reg1*norm(x, 2).^2 + reg2*norm(x, 1);
fzFunc.muhy        = 1.0e-1;
[optsol8, output8] = constrOptAsgardSolver(fzFunc, B, cb, options, y0);
y_opt8             = optsol8.x_opt;
fx_val8            = norm(B*y_opt8 - cb, 2) + 0.5*reg1*norm(y_opt8, 2)^2 + reg2*norm(y_opt8, 1);
diff_y8            = norm(y_opt8 - y_cvx, 2)/nrm_y_cvx;
time8              = toc(time8);

%% Call the PAPA solver for strongly convex case and with restart.
options.Algorithm  = 'CP';
time9              = tic;
fzFunc.usDefFunc   = @(x, y, varargin) norm(B*x - cb, 2) + 0.5*reg1*norm(x, 2).^2 + reg2*norm(x, 1);
fzFunc.muhy        = 1.0e-2;
[optsol9, output9] = constrOptAsgardSolver(fzFunc, B, cb, options, y0);
y_opt9             = optsol9.x_opt;
fx_val9            = norm(B*y_opt9 - cb, 2) + 0.5*reg1*norm(y_opt9, 2)^2 + reg2*norm(y_opt9, 1);
diff_y9            = norm(y_opt9 - y_cvx, 2)/nrm_y_cvx;
time9              = toc(time9);

%% Print the final results.
fprintf(['The objective values:       \n', ...
         '  +PAPA-SCVX(1.0)     = %3.20f\n  +PAPA-SCVX(0.1)     = %3.20f\n', ...
         '  +PAPA-SCVX(0.01)    = %3.20f\n  +PAPA-SCVX-RS(1.0)  = %3.20f\n', ...
         '  +PAPA-SCVX-RS(0.1)  = %3.20f\n  +PAPA-SCVX-RS(0.01) = %3.20f\n', ...
         '  +CP(1.0)            = %3.20f\n  +CP(0.1)            = %3.20f\n', ...
         '  +CP(0.01)           = %3.20f\n  +CVX                = %3.20f\n'], ...
         fx_val1, fx_val2, fx_val2, fx_val4, fx_val5, fx_val6, ...
         fx_val7, fx_val8, fx_val9, fx_cvx);
     
fprintf(['The solution differences:\n', ...
         '  +PAPA-SCVX(1.0)     = %3.20f\n  +PAPA-SCVX(0.1)     = %3.20f\n', ...
         '  +PAPA-SCVX(0.01)    = %3.20f\n  +PAPA-SCVX-RS(1.0)  = %3.20f\n', ...
         '  +PAPA-SCVX-RS(0.1)  = %3.20f\n  +PAPA-SCVX-RS(0.01) = %3.20f\n', ...
         '  +CP(1.0)            = %3.20f\n  +CP(0.1)            = %3.20f\n', ...
         '  +CP(0.01)           = %3.20f\n'], ...
         diff_y1, diff_y2, diff_y3, diff_y4, diff_y5, diff_y6, diff_y7, diff_y8, diff_y9);
     
fprintf(['The solution time in second:\n', ...
         '  +PAPA-SCVX(1.0)     = %3.20f\n  +PAPA-SCVX(0.1)     = %3.20f\n', ...
         '  +PAPA-SCVX(0.01)    = %3.20f\n  +PAPA-SCVX-RS(1.0)  = %3.20f\n', ...
         '  +PAPA-SCVX-RS(0.1)  = %3.20f\n  +PAPA-SCVX-RS(0.01) = %3.20f\n', ...
         '  +CP(1.0)            = %3.20f\n  +CP(0.1)            = %3.20f\n', ...
         '  +CP(0.01)           = %3.20f\n  +CVX                = %3.20f\n'], ...
         time1, time2, time3, time4, time5, time6, time7, time8, time9, time_cvx);

%% Plot the outputs.
fx_min  = min([fx_val1; fx_val2; fx_val3; fx_val4; fx_val5; fx_val6; fx_val7; fx_val8; fx_val9; fx_cvx]);
YMatrix{run} = [abs(output1.hist.usedef  - fx_min)/max(1, abs(fx_min)), ...
                abs(output2.hist.usedef  - fx_min)/max(1, abs(fx_min)), ...
                abs(output3.hist.usedef  - fx_min)/max(1, abs(fx_min)), ...
                abs(output4.hist.usedef  - fx_min)/max(1, abs(fx_min)), ...
                abs(output5.hist.usedef  - fx_min)/max(1, abs(fx_min)), ...
                abs(output6.hist.usedef  - fx_min)/max(1, abs(fx_min)), ...
                abs(output7.hist.usedef  - fx_min)/max(1, abs(fx_min)), ...
                abs(output8.hist.usedef  - fx_min)/max(1, abs(fx_min)), ...
                abs(output9.hist.usedef  - fx_min)/max(1, abs(fx_min))];

end
%%% End of two run

%% Generate a figure
createFigure1d(YMatrix{1}, YMatrix{2});

% PAPA v.1.0 by Quoc Tran-Dinh (quoctd@email.unc.edu)
% Copyright @ 2017 Department of Statistics and Operations Research (STOR)
%                The University of North Carolina at Chapel Hill (UNC)
% See the file LICENSE for full license information.

