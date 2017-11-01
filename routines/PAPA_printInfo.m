%% FUNCTION PAPA_printInfo()
%  PURPOSE: Print the header and information ...

% No output is required ...
if options.Verbosity <= 0, return; end
if options.Verbosity <= 1, fprintf('\n+ Iteration:     \n'); return; end

% Otherwise, print caption ...
fprintf('%s\n', repmat('*', 1, 66));
if strcmpi(options.Algorithm, 'PAPA')
    fprintf('       Proximal Alternating Penalized Algorithm \n');
	fprintf('                    **** PAPA **** \n');
elseif strcmpi(options.Algorithm, 'PAPA-RS')
    fprintf('     Restarting Proximal Alternating Penalized Algorithm \n');
	fprintf('                **** Restarting PAPA **** \n');
elseif strcmpi(options.Algorithm, 'PALPA')
    fprintf('     Proximal Alternating Linearized Penalized Algorithm \n');
	fprintf('                    **** PALPA **** \n');
elseif strcmpi(options.Algorithm, 'PALPA-RS')
    fprintf('Restarting Proximal Alternating Linearized Penalized Algorithm \n');
	fprintf('               **** Restarting PALPA **** \n');
elseif strcmpi(options.Algorithm, 'PAPA-SCVX')
    fprintf('     Proximal Alternating Linearized Penalized Algorithm \n');
    fprintf('           For Strongly Convex Objective Function \n');
	fprintf('                    **** PAPA-mu **** \n');
elseif strcmpi(options.Algorithm, 'PAPA-SCVX-RS')
    fprintf('Restarting Proximal Alternating Linearized Penalized Algorithm \n');
    fprintf('           For Strongly Convex Objective Function \n');
	fprintf('               **** Restarting PAPA-mu **** \n');
else
    fprintf('                     Unknown algorithm! \n');
end
fprintf('%s\n', repmat('*', 1, 66));
fprintf('PAPA v.1.0: A Proximal Alternating Penalty Optimization Solver.\n');
fprintf('       Copyright (c) 2017 Quoc Tran-Dinh (quoctd@email.unc.edu)\n');
fprintf('       Department of Statistics and Operations Research. \n');
fprintf('       The University of North Carolina at Chapel Hill (UNC). \n');
fprintf(['       See the file <a href="http://trandinhquoc.com">', ...
         'LICENSE</a> for full license information.\n']);

% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.