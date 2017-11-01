% FUNCTION: PAPA_saveHistory()
% PURPOSE:  Save the history w.r.t. the iterations if requires.
% 
% Level 1: Basic information: relative error, residual and objective values.
if options.saveHistMode > 0
    output.hist.rel_schg(  iter, 1) = rel_schg;
    output.hist.rel_pfeas( iter, 1) = rel_pfeas;
    output.hist.fx_val(    iter, 1) = fx_val;
    output.hist.time_it(   iter, 1) = toc(time_it);
end

% Level 2: Absolute error and residual.
if options.saveHistMode > 1
    output.hist.abs_schg(  iter, 1) = abs_schg;
    output.hist.abs_pfeas( iter, 1) = abs_pfeas;
end

% Level 3: Parameters.
if options.saveHistMode > 2
   output.hist.beta(iter, 1)  = beta;
   output.hist.gamma(iter, 1) = gamma;
   output.hist.tau( iter, 1)  = tau;
end
     
% Level 4: User define.
if options.saveHistMode > 3
    if isfield(objFunc, 'usDefFunc')
        output.hist.usedef(iter, :) = objFunc.usDefFunc(x_cur, y_cur);
    end
end

% Level 4: Save solution
if options.saveHistMode > 4
    if exist('x_cur'), output.hist.x_cur{iter} = x_cur; end
    if exist('y_cur'), output.hist.y_cur{iter} = y_cur; end
end

% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.