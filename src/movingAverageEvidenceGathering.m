function state = movingAverageEvidenceGathering(state, params)
% calculates a moving average of M as outcome of the EGP for the online
% stopping criteria framework

for j = 1:length(state.PI)
    M = state.(state.PI{j});
    if ~isfield(state, 'moving') || ...
            size(state.moving,1) < length(state.PI) || ...
            any(isnan(state.moving(j,[1 3])))
        % first generation: initialize moving average with first value and
        % initial guess of the noise
        if ~isnan(M(end))
            state.moving(j,:) = [M(end) NaN 1];
        else
            state.moving(j,:) = [NaN NaN NaN];
        end;
    elseif ~isnan(M(end))
        state.moving(j,:) = [mean([M(end) state.moving(j,1)]) ...
            std([M(end) state.moving(j,1)]) state.moving(j,3)+1];
    end;
end;