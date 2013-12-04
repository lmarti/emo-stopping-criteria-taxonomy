function state = rInd(PFt, state, params)
% implementation of the R indicator for the use in the 
% online stopping criteria framework
%
% ---
% Version: $Id: rInd.m 97 2010-10-09 16:33:46Z wagner $

% R indicator works in the objective space
currentFront = PFt.PF;

% check input
if ~isfield(state, 'preGen')
    % first generation: no values without reference generation
    state.r = NaN;
    return;
end

% compute R indicator values based on the current generation
n = length(state.preGen);
values = zeros(n,1);
for i = 1:n
    values(i) = rIndicator(state.preGen{i}, currentFront, ...
        params.idealPoint, params.refPoint);
end;
state.r = values;