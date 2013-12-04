function state = convergenceMetric(PFt, state, params)
% implementation of the convergence metric (also known as generational
% distance, GD) for the use in the online stopping criteria framework

% CM works in the objective space
currentFront = PFt.PF;

% compute cm for current generation
dist = zeros(size(currentFront,1),size(params.refSet,1));
for j = 1:size(dist,1)
    for k = 1:size(dist,2)
        dist(j,k) = norm(currentFront(j,:)-params.refSet(k,:));                
    end;
end;
cm = mean(min(dist, [], 2));

% update vector of stored indicator values
if isfield(state, 'cm')
    cm = cm ./ state.cmInit;
    state.cm = updateVector(cm, state.cm, params.tmem);
else
    state.cmInit = cm;
    state.cm = updateVector(1, [], params.tmem);    
end;