function state = consolidationRatio(PFt, state, params)
% implementation of the utility-based consolidation ratio for the use in 
% the online stopping criteria framework

% check input
if ~isfield(state, 'archive') || length(state.archive) < params.tmem
    % no values without archive of reference generation (t - tmem)
    state.cr = NaN;
    return;
end

% compute CR for current generation
ind = ismember(state.archive{1}, state.archive{params.tmem}, 'rows');
cr = sum(ind)./size(state.archive{params.tmem},1);

% update vector of stored indicator values
if ~isfield(state, 'crPreGen')
    % first generation: initialize crPreGen and iteration
    state.crPreGen = updateVector(cr, [], params);
    state.iteration = params.tmem;
    % no values without reference generation
    state.cr = NaN;
    state.iteration = state.iteration + 1;
    return;
else
    state.crPreGen = updateVector(cr, state.crPreGen, params.tmem);
    state.iteration = state.iteration + 1;
end

if cr > 0.5 && length(state.crPreGen) == params.tmem
    % compute average utility
    u = (cr - state.crPreGen(1))./params.tmem;
    if ~isfield(state, 'uInit')
        % adaptively set threshold
        idx = strmatch('cr', state.PI, 'exact');
        params.epsilon(:,idx) = (u./state.iteration)./params.F;
    end;
    % update vector of stored indicator values
    if any(isnan(state.cr))
        state.cr = updateVector(u, [], params.tmem);
    else
        state.cr = updateVector(u, state.cr, params.tmem);
    end;
else
    state.cr = NaN;
end;