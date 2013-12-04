function state = hypervolumeIndicator(PFt, state, params)
% implementation of the hypervolume indicator for the use in the
% online stopping criteria framework

% HV works in the objective space
currentFront = PFt.PF;

% check input
if ~isfield(state, 'hvPreGen')
    % first generation: initialize hvPreGen
    refVal = hv(currentFront', params.refPoint);
    state.hvPreGen = updateVector(refVal, [], params);
    % no values without reference generation
    state.hv = NaN;
    return;
end

% compute hv for current generation
refVal = hv(currentFront', params.refPoint);

% compare indicator values of preceding generations to the current generation
state.hv = refVal - state.hvPreGen;

% update vector of stored indicator values
state.hvPreGen = updateVector(refVal, state.hvPreGen, params.tmem);