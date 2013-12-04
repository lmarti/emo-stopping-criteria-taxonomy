function state = directTransfer(state, params)
% directly transfers the value M to the outcome of the EGP for the online
% stopping criteria framework

for j = 1:length(state.PI)
    M = state.(state.PI{j});    
    state.direct(j,:) = [M(end) NaN 1];    
end;