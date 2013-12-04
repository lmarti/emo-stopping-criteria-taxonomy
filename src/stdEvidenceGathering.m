function state = stdEvidenceGathering(state, params)
% computing the standard deviation of the values in M as outcome of the 
% EGP for the online stopping criteria framework

for j = 1:length(state.PI)
    M = state.(state.PI{j});
    if length(M) < params.tmem || any(isnan(M))       
        state.std(j,:) = [NaN NaN NaN];
    else
        val = std(M);
        btrpResults = bootstrp(params.btrpSamples,@std,M);
        error = sqrt(sum((val-btrpResults).^2.)/(params.btrpSamples-1));
        state.std(j,:) = [val error params.tmem-1];
    end;       
end;