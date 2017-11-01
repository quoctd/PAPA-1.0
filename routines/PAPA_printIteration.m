% FUNCTION: PAPA_printIteration()
% PURPOSE:  Print the iteration if required.
% 
% Not output is printed.
if options.Verbosity <= 1, return; end

% Print the header.
if mod(iter, 10*options.PrintStep) == 1 || iter == 1
    fprintf('%s\n', repmat('-', 1, options.PrintLength)); 
    fprintf([' Iter| RelPGap | RelSchg |', ...
             '  Beta  | Gamma  |   Tau  |    F(x)\n']);
    fprintf('%s\n', repmat('-', 1, options.PrintLength));
end

% Print the values: iter -> rel_pfeas -> rel_dfeas -> rel_schg 
%                   beta -> gamma -> tau -> fx_val
if mod(iter, options.PrintStep) == 0 || iter == 1
    fprintf('%5d| %3.2e| %3.2e| %3.1e| %3.1e| %3.1e| %3.5e \n', ...
             iter, rel_pfeas, rel_schg, beta, gamma, tau, fx_val);
end

% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.