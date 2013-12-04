function state = epsInd(PFt, state, params)
% implementation of the epsilon indicator for the use in the 
% online stopping criteria framework
%
% ---
% Version: $Id: epsInd.m 97 2010-10-09 16:33:46Z wagner $

% Epsilon indicator works in the objective space
currentFront = PFt.PF;

% check input
if ~isfield(state, 'preGen')
    % first generation: no values without reference generation
    state.epsilon = NaN;
    return;
end

% compute epsilon indicator values based on the current generation
n = length(state.preGen);
values = zeros(n,1);
for i = 1:n
    values(i) = epsilonIndicator(state.preGen{i}, currentFront);
end;
state.epsilon = values;