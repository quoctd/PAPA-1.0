% The problem is
%                   min  g(x) + r(y)
%                   s.t. B*y - x = 0.
% where g(x) = lambda*|x|_1 and r(y) = 0.5*|A*y - b|^2.     
%
% To run this, we need to download UNLOCBOX for computing the proximal 
% operator of the TV-norm. See https://epfl-lts2.github.io/unlocbox-html/

% This controls which algorithms to run (1 = run, or 0 = not run).
isRunAlg                = [1, 1, 1, 1, 1, 1, 1, 1];

file_name               = 'im1_MRI_hip';
lambda                  = 4.0912e-04; % Without noise and with noise (HIP)

data_opt                = 4;
imaddress               = [file_name, '.jpg']; %
small_data_size         = [16, 16];   % optional, only required if 3 is chosen above
nresize                 = [256, 256]; % Rescale to work with small image.
[x_original, ind, m, n] = get_dataset(data_opt, small_data_size, imaddress, nresize);                                
x_original              = real(x_original);
 
% Define the A=subsampled fourier and D=difference operators
Aoper                 = @(X) fft2fwd(X, ind, m, n);
AToper                = @(X) fft2adj_rectangular_new(X, ind, m, n);
Boper                 = @(X) gradOperator( X );
BToper                = @(X) divOperator( X, n );

% The observed data and regularization parameter.
cb                    = Aoper(x_original);

% Define an initial input.
Y0                    = zeros(size(x_original));
X0                    = Boper(Y0);

% Compute the Lipschitz constant of ry and the norm of B.
LipsR                 = 1.0;
LB_bar                = 8.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set the optional parameters.
options                = PAPA_OptimSet([]);
options.PrintStep      = 30;
options.isStoppingCond = 0;
options.MaxIters       = 200;
options.nRestart       = 100;
options.PwMaxIters     = 30;
options.saveHistMode   = 4;
options.Verbosity      = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            CALL SOLVERS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************************************************************
%% Define the objective function and its proximal operators.
%**************************************************************************
px                    = size(X0);
py                    = size(Y0);
objFunc.px            = px;
objFunc.py            = py;
objFunc.nc            = px;

soft_threshold_oper   = @(X, gamma) sign(X).*max(abs(X) - gamma, 0);
objFunc.gxProxOper    = @(X, gamma, varargin) ( soft_threshold_oper(X, lambda*gamma) );
objFunc.hyProxOper    = @(Y, gamma, varargin) ( Y );
objFunc.ryGradOper    = @(Y, varargin) ( AToper(Aoper(Y) - cb) );

objFunc.gxFunc        = @(X, varargin) ( lambda*norm(X(:), 1) );
objFunc.hyFunc        = @(Y, varargin) ( 0 );
objFunc.ryFunc        = @(Y, varargin) ( 0.5*norm( Aoper(Y) - cb, 'fro').^2 );
objFunc.ryLips        = LipsR;

linConstr.Aoper       = @(X, varargin) ( -X );
linConstr.AToper      = @(X, varargin) ( -X );
linConstr.Boper       = @(Y, varargin) ( Boper( Y ) );
linConstr.BToper      = @(Y, varargin) ( BToper(Y) );
linConstr.dlStarProx  = @(U, gamma, varargin) ( U );
linConstr.dlStarFunc  = @(U, varargin) ( norm(U, 'fro') );

% User define function to compute the objective value.
objFunc.usDefFunc     = @(X, Y, varargin) 0.5*norm( Aoper(Y) - cb, 'fro').^2 + lambda*norm(vec(Boper(Y)), 1);

% Set other parameters.
LA_bar                = 1;
linConstr.LA_bar      = LA_bar;
linConstr.LB_bar      = LB_bar;

%% Define the parameters for PAPA.
muhy                  = 0.5;
beta1                 = 2.0*sqrt(LB_bar);
TV_prox_Iter1         = 25;
TV_prox_Iter2         = 100;

%**************************************************************************
%% Call our PAPA solver - non-strong convexity and no-restart.
%**************************************************************************
if isRunAlg(1)
    options.Algorithm     = 'PAPA';
    options.isRestart     = 0;

    linConstr.beta1       = beta1;
    time1                 = tic;
    [optsol1, output1]    = nscvxPapa3Solver(objFunc, linConstr, X0, Y0, options);
    time1                 = toc(time1);
    % Get the output.
    Xopt1                 = abs(optsol1.y_opt);
    fx_val1               = 0.5*norm( Aoper(Xopt1) - cb, 'fro').^2 + lambda*norm(vec(Boper(Xopt1)), 1);
    sol_diff1             = norm(x_original - Xopt1, 'fro')/norm(x_original, 'fro');
    psnr1                 = psnr(Xopt1, x_original);
else
    fx_val1 = nan; sol_diff1 = nan; psnr1 = nan; time1 = nan;
end
%**************************************************************************
%% Call our PAPA solver - non-strong convexity and with restart.
%**************************************************************************
if isRunAlg(2)
    options.Algorithm     = 'PAPA-RS';
    options.isRestart     = 1;

    linConstr.beta1       = beta1;
    time2                 = tic;
    [optsol2, output2]    = nscvxPapa3Solver(objFunc, linConstr, X0, Y0, options);
    time2                 = toc(time2);
    % Get the output.
    Xopt2                 = abs(optsol2.y_opt);
    fx_val2               = 0.5*norm( Aoper(Xopt2) - cb, 'fro').^2 + lambda*norm(vec(Boper(Xopt2)), 1);
    sol_diff2             = norm(x_original - Xopt2, 'fro')/norm(x_original, 'fro');
    psnr2                 = psnr(Xopt2, x_original);
else
    fx_val2 = nan; sol_diff2 = nan; psnr2 = nan; time2 = nan; Xopt2 = [];
end

%**************************************************************************
%% Call our PAPA solver - strong convexity and no-restart.
%**************************************************************************
if isRunAlg(3)
    options.Algorithm      = 'PAPA-SCVX';
    options.isRestart      = 0;

    linConstr.beta1        = 2*LB_bar/muhy;
    time3                  = tic;
    [optsol3, output3]     = scvxPapa3Solver(objFunc, linConstr, X0, Y0, options);
    time3                  = toc(time3);
    % Get the output.
    Xopt3                  = abs(optsol3.y_opt);
    fx_val3                = 0.5*norm( Aoper(Xopt3) - cb, 'fro').^2 + lambda*norm(vec(Boper(Xopt3)), 1);
    sol_diff3              = norm(x_original - Xopt3, 'fro')/norm(x_original, 'fro');
    psnr3                  = psnr(Xopt3, x_original);
else
    fx_val3 = nan; sol_diff3 = nan; psnr3 = nan; time3 = nan; Xopt3 = [];
end


%**************************************************************************
%% Call our PAPA solver - strong convexity and with restart.
%**************************************************************************
if isRunAlg(4)
    options.Algorithm      = 'PAPA-SCVX-RS';
    options.isRestart      = 1;

    linConstr.beta1        = 2*LB_bar/muhy;
    time4                  = tic;
    [optsol4, output4]     = scvxPapa3Solver(objFunc, linConstr, X0, Y0, options);
    time4                  = toc(time4);
    % Get the output.
    Xopt4                  = abs(optsol4.y_opt);
    fx_val4                = 0.5*norm( Aoper(Xopt4) - cb, 'fro').^2 + lambda*norm(vec(Boper(Xopt4)), 1);
    sol_diff4              = norm(x_original - Xopt4, 'fro')/norm(x_original, 'fro');
    psnr4                  = psnr(Xopt4, x_original);
else
    fx_val4 = nan; sol_diff4 = nan; psnr4 = nan; time4 = nan; Xopt4 = [];
end
%**************************************************************************
%% Call Vu-Condat's algorithm without tuning.
%**************************************************************************
objFunc2.nx           = size(Y0);
soft_threshold_oper   = @(X, gamma) sign(X).*max(abs(X) - gamma, 0);
objFunc2.gxProxOper   = @(X, gamma, varargin) ( soft_threshold_oper(X, lambda*gamma) );
objFunc2.hxProxOper   = @(X, gamma, varargin) ( X );
objFunc2.rxGradOper   = @(X, varargin) ( AToper(Aoper(X) - cb) );

objFunc2.gxFunc       = @(X, varargin) ( lambda*norm(X(:), 1) );
objFunc2.hxFunc       = @(X, varargin) ( 0 );
objFunc2.rxFunc       = @(X, varargin) ( 0.5*norm( Aoper(X) - cb, 'fro').^2 );
objFunc2.rxLips       = LipsR;

linOper2.Boper        = @(Y, varargin) ( Boper( Y ) );
linOper2.BToper       = @(Y, varargin) ( BToper(Y) );
% User define function to compute the objective value.
objFunc2.usDefFunc    = @(X, Y, varargin) 0.5*norm( Aoper(X) - cb, 'fro').^2 + lambda*norm(vec(Boper(X)), 1);

% Set other parameters.
linOper2.LB_bar       = LB_bar;

%% Call Vu-Condat's algorithm with tuned parameter.
%**************************************************************************
if isRunAlg(5)
    param.tau             = 0.089/LipsR;
    param.sigma           = (1/param.tau - 0.5*LipsR)/LB_bar;
    param.theta           = 1.0;
    options.Algorithm     = 'VU-CONDAT';
    time6                 = tic;
    [optsol6, output6]    = nscvxVuCondatSolver(objFunc2, linOper2, Y0, options, param);
    time6                 = toc(time6);
    Xopt6                 = abs(optsol6.x_opt);
    fx_val6               = 0.5*norm( Aoper(Xopt6) - cb, 'fro').^2 + lambda*norm(vec(Boper(Xopt6)), 1);
    sol_diff6             = norm(x_original - Xopt6, 'fro')/norm(x_original, 'fro');
    psnr6                 = psnr(Xopt6, x_original);
else
    fx_val6 = nan; sol_diff6 = nan; psnr6 = nan; time6 = nan; Xopt6 = [];
end
%**************************************************************************
%% Call Accelerated Proximal-Gradient algorithm without restart.
%**************************************************************************
% Define the proximal operator.
param_tv.verbose = 0;
param_tv.tol     = 1e-10;
param_tv.maxit   = TV_prox_Iter1; 
TvNormProx       = @(X, gamma) prox_tv(X, gamma*lambda, param_tv);

% Get the problem input data.
objFunc3.nx           = size(Y0);
objFunc3.fxGradOper   = @(X, varargin) ( AToper( Aoper(X) - cb ) );
objFunc3.gxProxOper   = @(X, gamma, varargin) ( TvNormProx(X, gamma) );
objFunc3.fxFunc       = @(X, varargin) 0.5*norm( Aoper(X) - cb, 'fro').^2;
objFunc3.gxFunc       = @(X, varargin) lambda*norm(vec(Boper(X)), 1);
objFunc3.fxLips       = LipsR;
% User define function to compute the objective value.
objFunc3.usDefFunc    = @(X, Y, varargin) 0.5*norm( Aoper(X) - cb, 'fro').^2 + lambda*norm(vec(Boper(X)), 1);

% Call the solver.
if isRunAlg(6)
    options.isRestart     = 0;
    options.Algorithm     = 'ACC-PROX-GRAD';
    time7                 = tic;
    [optsol7, output7]    = AccProxGradSolver(objFunc3, Y0, options);
    time7                 = toc(time7);
    Xopt7                 = abs(optsol7.x_opt);
    fx_val7               = 0.5*norm( Aoper(Xopt7) - cb, 'fro').^2 + lambda*norm(vec(Boper(Xopt7)), 1);
    sol_diff7             = norm(x_original - Xopt7, 'fro')/norm(x_original, 'fro');
    psnr7                 = psnr(Xopt7, x_original);
else
    fx_val7 = nan; sol_diff7 = nan; psnr7 = nan; time7 = nan; Xopt7 = [];
end
%**************************************************************************
%% Call Accelerated Proximal-Gradient algorithm with restart.
%**************************************************************************
% Call the solver.
if isRunAlg(7)
    options.isRestart     = 1;
    options.Algorithm     = 'ACC-PROX-GRAD-RS';
    time8                 = tic;
    [optsol8, output8]    = AccProxGradSolver(objFunc3, Y0, options);
    time8                 = toc(time8);
    Xopt8                 = abs(optsol8.x_opt);
    fx_val8               = 0.5*norm( Aoper(Xopt8) - cb, 'fro').^2 + lambda*norm(vec(Boper(Xopt8)), 1);
    sol_diff8             = norm(x_original - Xopt8, 'fro')/norm(x_original, 'fro');
    psnr8                 = psnr(Xopt8, x_original);
else
    fx_val8 = nan; sol_diff8 = nan; psnr8 = nan; time8 = nan; Xopt8 = [];
end
%**************************************************************************
%% Call Accelerated Proximal-Gradient algorithm with restart.
%**************************************************************************
% Define the proximal operator.
param_tv.maxit   = TV_prox_Iter2; 
TvNormProx       = @(X, gamma) prox_tv(X, gamma*lambda, param_tv);
%TvNormProx = @(X, gamma) chambolle_prox_TV_stop(X, gamma*lambda, TV_prox_Iter2, 1e-10);
objFunc3.gxProxOper   = @(X, gamma, varargin) ( TvNormProx(X, gamma) );

% Call the solver.
if isRunAlg(8)
    options.isRestart     = 1;
    options.Algorithm     = 'ACC-PROX-GRAD2';
    time9                 = tic;
    [optsol9, output9]    = AccProxGradSolver(objFunc3, Y0, options);
    time9                 = toc(time9);
    Xopt9                 = abs(optsol9.x_opt);
    fx_val9               = 0.5*norm( Aoper(Xopt9) - cb, 'fro').^2 + lambda*norm(vec(Boper(Xopt9)), 1);
    sol_diff9             = norm(x_original - Xopt9, 'fro')/norm(x_original, 'fro');
    psnr9                 = psnr(Xopt9, x_original);
else
    fx_val9 = nan; sol_diff9 = nan; psnr9 = nan; time9 = nan; Xopt9 = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           PRINT OUTPUTS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fprintf(' -- The objective value: \n');
fprintf(['     +PAPA-nscvx     = %3.10f\n     +PAPA-nscvx-rs  = %3.10f\n', ...
         '     +PAPA-scvx      = %3.10f\n     +PAPA-scvx-rs   = %3.10f\n', ...
         '     +Vu-Condat-tune = %3.10f\n     +AccProxGrad-25 = %3.10f\n', ...
         '     +AccProxGrad-rs = %3.10f\n     +AccProxGrad-50 = %3.10f\n'], ...
         fx_val1, fx_val2, fx_val3, fx_val4, fx_val6, fx_val7, fx_val8, fx_val9);
fprintf(' -- The solution difference: \n');
fprintf(['     +PAPA-nscvx     = %3.10f\n     +PAPA-nscvx-rs  = %3.10f\n', ...
         '     +PAPA-scvx      = %3.10f\n     +PAPA-scvx-rs   = %3.10f\n', ...
         '     +Vu-Condat-tune = %3.10f\n     +AccProxGrad-25 = %3.10f\n', ...
         '     +AccProxGrad-rs = %3.10f\n     +AccProxGrad-50 = %3.10f\n'], ...
         sol_diff1, sol_diff2, sol_diff3, sol_diff4, ...
         sol_diff6, sol_diff7, sol_diff8, sol_diff9);
fprintf(' -- The PSNR: \n');
fprintf(['     +PAPA-nscvx     = %3.10f\n     +PAPA-nscvx-rs  = %3.10f\n', ...
         '     +PAPA-scvx      = %3.10f\n     +PAPA-scvx-rs   = %3.10f\n', ...
         '     +Vu-Condat-tune = %3.10f\n     +AccProxGrad-25 = %3.10f\n', ...
         '     +AccProxGrad-rs = %3.10f\n     +AccProxGrad-50 = %3.10f\n'], ...
         psnr1, psnr2, psnr3, psnr4, psnr6, psnr7, psnr8, psnr9);
fprintf(' -- The solution time: \n');
fprintf(['     +PAPA-nscvx     = %3.4f\n     +PAPA-nscvx-rs  = %3.4f\n', ...
         '     +PAPA-scvx      = %3.4f\n     +PAPA-scvx-rs   = %3.4f\n', ...
         '     +Vu-Condat-tune = %3.4f\n     +AccProxGrad-25 = %3.4f\n', ...
         '     +AccProxGrad-rs = %3.4f\n     +AccProxGrad-50 = %3.4f\n'], ...
         time1, time2, time3, time4, time6, time7, time8, time9);
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               PLOTTING                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the results.
styles = {'-', '--', '-.', ':'};
colors = {'m', 'r', 'b', 'g', 'k', [0.5, 0.5, 0], [0.5, 0, 1], [0.25, 0.75 0], [0, 0.25, 0.75]};
myplot = @(x, s, c) semilogy(x, 'MarkerSize', 3, 'Color', colors{c}, 'LineWidth', 2, 'LineStyle', styles{s});
myabs  = @(x) (x);
fx_min = 0;

if options.saveHistMode > 3
    figure(1);
    title('The real objective residual');
    strlegend = {};
    if isRunAlg(1), myplot(myabs(output1.hist.usedef - fx_min)/max(1, abs(fx_min)), 1, 1); hold on; strlegend = {strlegend{:}, 'PAPA'}; end
    if isRunAlg(2), myplot(myabs(output2.hist.usedef - fx_min)/max(1, abs(fx_min)), 2, 2); hold on; strlegend = {strlegend{:}, 'PAPA-rs'}; end
    if isRunAlg(3), myplot(myabs(output3.hist.usedef - fx_min)/max(1, abs(fx_min)), 3, 3); hold on; strlegend = {strlegend{:}, 'scvx-PAPA'}; end
    if isRunAlg(4), myplot(myabs(output4.hist.usedef - fx_min)/max(1, abs(fx_min)), 4, 4); hold on; strlegend = {strlegend{:}, 'scvx-PAPA-rs'}; end
    if isRunAlg(5), myplot(myabs(output6.hist.usedef - fx_min)/max(1, abs(fx_min)), 2, 5); hold on; strlegend = {strlegend{:}, 'Vu-Condat-tuned'}; end
    if isRunAlg(6), myplot(myabs(output7.hist.usedef - fx_min)/max(1, abs(fx_min)), 3, 6); hold on; strlegend = {strlegend{:}, 'AcProxGrad-25'}; end
    if isRunAlg(7), myplot(myabs(output8.hist.fx_val - fx_min)/max(1, abs(fx_min)), 4, 7); hold on; strlegend = {strlegend{:}, 'AcProxGrad-rs'}; end
    if isRunAlg(8), myplot(myabs(output9.hist.fx_val - fx_min)/max(1, abs(fx_min)), 1, 8); hold on; strlegend = {strlegend{:}, 'AcProxGrad-50'}; end
    legend(strlegend);
end

%% Plot the image.
figure(2);
colormap('gray');
subplot(3,3,1); imagesc(x_original); hold on; xlabel('Original');
if isRunAlg(1), subplot(3,3,2); imagesc(Xopt1); hold on; xlabel('PAPA');            end
if isRunAlg(2), subplot(3,3,3); imagesc(Xopt2); hold on; xlabel('PAPA-rs');         end
if isRunAlg(3), subplot(3,3,4); imagesc(Xopt3); hold on; xlabel('scvx-PAPA');       end
if isRunAlg(4), subplot(3,3,5); imagesc(Xopt4); hold on; xlabel('scvx-PAPA-rs');    end    
if isRunAlg(5), subplot(3,3,6); imagesc(Xopt6); hold on; xlabel('Vu-Condat-tuned'); end
if isRunAlg(6), subplot(3,3,7); imagesc(Xopt7); hold on; xlabel('AcProxGrad-25');      end
if isRunAlg(7), subplot(3,3,8); imagesc(Xopt8); hold on; xlabel('AcProxGrad-rs');   end
if isRunAlg(8), subplot(3,3,9); imagesc(Xopt9); hold on; xlabel('AcProxGrad-50');    end

% PAPA v.1.0 by Quoc Tran-Dinh (quoctd@email.unc.edu)
% Copyright @ 2017 Department of Statistics and Operations Research (STOR)
%                The University of North Carolina at Chapel Hill (UNC)
% See the file LICENSE for full license information.