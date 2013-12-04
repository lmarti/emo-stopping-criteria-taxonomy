function state = maximumCrowdingIndicator(PFt, state, params)
% implementation of the maximum crowding distance for the use in the 
% online stopping criteria framework

% maxCD works in the objective space
currentFront = PFt.PF;

% compute maxCD for current generation
dist = zeros(size(currentFront,1),size(currentFront,2));
preI = max(0:size(currentFront,1)-1,1);
postI = min(2:size(currentFront,1)+1,size(currentFront,1));
for j = 1:size(currentFront, 2)
    [val idx] = sort(currentFront(:,j));
    dist(idx,j) = val(postI) - val(preI);    
end;
maxCD = max(sum(dist,2));

% update vector of stored indicator values
if isfield(state, 'maxCD')
    state.maxCD = updateVector(maxCD, state.maxCD, params.tmem);
else
    state.maxCD = updateVector(maxCD, [], params.tmem);
end;