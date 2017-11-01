% FUNCTION: param = PAPA_OptimSet(varargin)
%
% PURPOSE:  Generate optional parameters for algorithms.
% 
% CALL: 
%    help    PAPA_OptimSet: see optional parameters.   
%    param = PAPA_OptimSet() to see the default values.
%    param = PAPA_OptimSet('Name', Value, ...) to set new param.
%    param = PAPA_OptimSet(param, 'Name', Value, ...) to add new param.
%
% NOTES:
%    MaxIters        : The maximum number of iterations.
%    RelTolX         : The relative tolerance for the search direction.
%    RelTolFun       : The relative tolerance for the objective value.
%    RelTolFeas      : The relative tolerance for the feasibility gap.
%    TerminationType : Type of termination conditions.', char(13), ...
%                      (0 = default, 1 = norm of the search
%                      direction, 2 - objective value change, 3 = 1 + 2).
%    isFxEval        : Evaluating the objective value if required.
%    Verbosity       : Output printing level.
%    isRestart       : Perform restarting step.
%    nRestart        : The number of iterations to perform Restarting.
%    saveHistMode    : Save the history information (0 = non and 1 = details, 
%                      2 = include the iterative vector).
%    PrintStep       : Number of iterations for printing output (50).
%    PwMaxIters      : Maximum number of iterations for Power method.
%    PwRelTol        : The relative accuracy for Power method.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Department of Statistics and Operations Research
%    The University of North Carolina at Chapel Hill (UNC).
%    Date: 09.20.2017.
%    Last modified: 10.15.2017.
%    Contact: quoctd@email.unc.edu
%
function param = PAPA_OptimSet(varargin)

% Define the default param.
defaultopt = struct( ...
                     'Algorithm',       'PAPA',   ...
                     'MaxIters',        5000,       ...
                     'RelTolX',         1e-5,       ...
                     'RelTolFun',       1e-4,       ...
                     'RelTolFeas',      1e-5,       ...
                     'isFxEval',        true,       ...
                     'isRestart',       true,       ...
                     'nRestart',        50,        ...
                     'TerminationType', 1,          ...
                     'Verbosity',       2,          ...
                     'saveHistMode',    1,          ...
                     'PrintStep',       200,         ...
                     'PrintLength',     65,         ...
                     'TradeOffFact',    1.0,        ...
                     'lbScvxParam',     1e-5,       ...
                     'isStoppingCond',  1,          ...
                     'PwMaxIters',      30,         ...
                     'PwRelTol',        1e-4        ...
                   );
             
% In case no input or empty input.
if nargin < 1,        
    param = defaultopt; 
    printOptions();
    return; 
end
if isempty(varargin{1})
    param = defaultopt; 
    return; 
end

% In case varargin{1} is a structure.
if isstruct(varargin{1})
    param = varargin{1};
    if nargin > 1
        inputs = varargin(2:end);
    else
        inputs = [];
    end
else
   param = defaultopt; 
   inputs = varargin;
end
if mod(length(inputs), 2) ~= 0
    error('Optional parameters must be in pair (Name, Value)!');
end

% Now, add the new param.
param = addNewOptions(param, inputs);

% Check if inputs are correct.
param = checkfield(param, defaultopt, 'Algorithm' );
param = checkfield(param, defaultopt, 'MaxIters' );
param = checkfield(param, defaultopt, 'RelTolX' );
param = checkfield(param, defaultopt, 'RelTolFun' );
param = checkfield(param, defaultopt, 'RelTolFeas' );
param = checkfield(param, defaultopt, 'isFxEval' );
param = checkfield(param, defaultopt, 'Verbosity' );
param = checkfield(param, defaultopt, 'TerminationType');
param = checkfield(param, defaultopt, 'saveHistMode');
param = checkfield(param, defaultopt, 'PrintStep' );
param = checkfield(param, defaultopt, 'PrintLength' );
param = checkfield(param, defaultopt, 'PwMaxIters' );
param = checkfield(param, defaultopt, 'PwRelTol' );
param = checkfield(param, defaultopt, 'isRestart' );
param = checkfield(param, defaultopt, 'nRestart' );
param = checkfield(param, defaultopt, 'TradeOffFact' );
param = checkfield(param, defaultopt, 'lbScvxParam' );
param = checkfield(param, defaultopt, 'isStoppingCond' );

% FUNCTION: outopt = checkfield(outopt, default, fieldname) 
% PURPOSE:  This function checks the field of optional parameters.
function outopt = checkfield(outopt, default, fieldname)
    
    fieldvalue = getfield( default, fieldname );
    if isfield(outopt, fieldname )
        optfieldval = getfield(outopt, fieldname );
        if isempty( optfieldval )
            outopt = setfield( outopt, fieldname, fieldvalue );
        end;
    else
        outopt = setfield( outopt, fieldname, fieldvalue );
    end;

% FUNCTION: param = addNewOptions(param, inputs)
%PURPOSE:  Add new param to the structure.
function param = addNewOptions(defin, inputs)

param = defin;
if isempty(inputs), return; end

for k=1:length(inputs)/2
    name  = inputs{2*k-1};
    value = inputs{2*k};
    if ~ischar(name)
        error('Name of an option must be a string!');
    end
    if ~isfield(param, name)
        error('This optional parameter does not exist! Add a new one.');
    end
    param = setfield(param, name, value);
end

% FUNCTION: printOptions()
% PURPOSE:  Print the param.
function printOptions()

fprintf('  MaxIters      : The maximum number of iterations.\n');
fprintf(['  RelTolX       : The relative tolerance for the search ', ...
         'direction.\n']);
fprintf(['  RelTolFun     : The relative tolerance for the objective ', ...
         'values.\n']);
fprintf('  RelTolFeas    : The relative tolerance for the feasibility gap.\n');
fprintf('  isFxEval      : Evaluating the objective value if required.\n');
fprintf( '  Verbosity     : Option for printing output.\n');
fprintf( '  isRestart     : Perform a Restarting step.\n');
fprintf( '  nRestart      : Number of itertations to perform restarting.\n');
fprintf(['  TerminateType : Type of termination conditions.', char(13), ...
         '                  (0 = default, 1 = norm of the search', ...
         ' direction, 2 - objective value change, 3 = 1 + 2).\n']);
fprintf(['  saveHistMode  : Save the history information mode (0 = non', ...
         ', 1 = details, 2 - iterative vector).\n']);
fprintf('  PrintStep     : Number of iterations for printing output (50).\n');
fprintf('  PwMaxIters    : Maximum number of iterations for Power method.\n');
fprintf('  PwRelTol      : The relative accuracy for Power method.\n');
     
% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.