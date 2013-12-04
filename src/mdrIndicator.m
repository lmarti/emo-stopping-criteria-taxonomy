function state = mdrIndicator(nonDomCurrent, state, params)
%MDRINDICATOR Calculates the mutual domination rate indicator.
%   
%   Parameters:
%    * nonDomCurrent, non-dominated front and set as a matrix of n rows and M+D columns.
%    * state, state of the OSC
%    * params, struct providing the size of the time window (params.tmem)
%
%   Return:
%   * state.mdr, values of the indicator.
%
%   See:
%   Marti, L., Garcia, J., Berlanga, A., & Molina, J. M. (2009). 
%   An Approach to Stopping Criteria for Multi-Objective Optimization 
%   Evolutionary Algorithms: The MGBM Criterion. In 2009 IEEE Conference
%   on Evolutionary Computation (CEC 2009), Piscataway, New Jersey.
%
% ---
% Version: $Id: mdrIndicator.m 97 2010-10-09 16:33:46Z wagner $

% MDR works in the objective space
PFt = nonDomCurrent.PF;

% check input
if ~isfield(state, 'preGen')
    % first generation: no values without reference generation
    state.mdr = NaN;
    return;
end

% compute MDR
n = length(state.preGen);
mdr = zeros(n,1);
for k = 1:n 
    [normA normB] = normDelta(state.preGen{k}, PFt);   
    mdr(k) = normA - normB;  
end
state.mdr = mdr;