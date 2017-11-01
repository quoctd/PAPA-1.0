%% FUNCTION: PAPA_finalization
%  PURPOSE:  Finalizes the history vectors.

% Delete the rest of hist vectors if required.
% % Level 1: Basic information: relative error, residual and objective values.
% if options.saveHistMode > 0 && iter < options.MaxIters
%     output.hist.rel_schg(  iter+1:end, 1) = [];
%     output.hist.rel_pfeas( iter+1:end, 1) = [];
%     output.hist.fx_val(    iter+1:end, 1) = [];
% end

% % Level 2: Absolute error and residual.
% if options.saveHistMode > 1 && iter < options.MaxIters
%     output.hist.abs_schg(  iter+1:end, 1) = [];
%     output.hist.abs_pfeas( iter+1:end, 1) = [];
% end

% % Level 3: Parameters.
% if options.saveHistMode > 2 && iter < options.MaxIters
%    output.hist.eta(iter+1:end, 1)   = [];
%    output.hist.rho(iter+1:end, 1)   = [];
%    output.hist.beta(iter+1:end, 1)  = [];
%    output.hist.gamma(iter+1:end, 1) = [];
%    output.hist.tau( iter+1:end, 1)  = [];
% end

% Print the last iteration.
if options.Verbosity > 1
    % Print the values: iter -> rel_pfeas -> rel_dfeas -> rel_schg 
    %                   beta -> gamma -> tau -> fx_val
    if (mod(iter, options.PrintStep) == 0 || iter == 1 ) && iter < options.MaxIters
        fprintf('%5d| %3.2e| %3.2e| %3.1e| %3.1e| %3.1e| %3.5e \n', ...
             iter, rel_pfeas, rel_schg, beta, gamma, tau, fx_val);
    end
    fprintf('%s\n', repmat('*', 1, 66));
    fprintf('  End of algorithm. \n');
    fprintf('%s\n', repmat('*', 1, 66));
end

% If exceed the number of iteration.
if iter >= options.MaxIters
    output.status = 'Have not yet converged';
    output.msg    = 'Exceed the maximum number of iterations';
end

% Get other information.
output.iter      = iter;
output.fx_val    = fx_val;
output.rel_schg  = rel_schg;
output.rel_pfeas = rel_pfeas;
output.time      = toc(time1);

% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.