function state = updatePreGen(nonDomCurrent, state, params)
% update previous PF*

if ~isfield(state, 'preGen') || isempty(state.preGen)
    % initialize preGen in the beginning of the evolution
    state.preGen{1} = nonDomCurrent;
else
    n = length(state.preGen);
    if n < params.tmem
        state.preGen{n+1} = nonDomCurrent;
    else
        for k = 1:n-1;
            state.preGen{k} = state.preGen{k+1};
        end;
        state.preGen{n} = nonDomCurrent;
    end
end;
