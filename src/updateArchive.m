function state = updateArchive(nonDomCurrent, state, params)
% update archives of non-dominated solutions

if ~isfield(state, 'archive') || isempty(state.archive)
    % initialize archive in the beginning of the evolution
    state.archive{1} = nonDomCurrent;
else
    n = length(state.archive);
    archive = [state.archive{n}; nonDomCurrent];
    archive = archive(paretofront(archive),:);
    if n < params.tmem        
        state.archive{n+1} = archive;
    else
        for k = 1:n-1;
            state.archive{k} = state.archive{k+1};
        end;
        state.archive{n} = archive;
    end
end;