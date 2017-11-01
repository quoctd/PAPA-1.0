%% FUNCTION: PAPA_initHistory()
%  PURPOSE:  Initialize the history of the iterations.

%%% Initialize some undefined paramters ...
fx_val = nan; rel_schg = nan; rel_pfeas = nan; gamma = 0; beta = 0; tau = 0;
output.msg = ''; output.status = ''; gy_val = nan; dl_val = 0;
x_cur  = []; y_cur = [];
    
% PAPA v.1.0: @copyright by Quoc Tran-Dinh, 2017.
%    Department of Statistics and Operations Research
%    University of North Carolina at Chapel Hill, NC, USA.
% See the file LICENSE for full license information.