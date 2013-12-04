params.objFunc = load('UF2_Ref.txt');
params.refPoint   = [5 5];
params.idealPoint = [0 0];
params.objFunc = 'UF2';
% Rudenko and Schoenauer - Stability measure
[stopGeneration state] = oscTaxonomy('nsga2', 'uf2', 2, 1000, 1, ...
    {'maxCD'}, {'std'}, {'threshold'}, params)
% Marti et al. - MGBM criterion
[stopGeneration state] = oscTaxonomy('nsga2', 'uf2', 2, 1000, 1, ...
    {'mdr'}, {'kalman'}, {'ciNormal'}, params)
% Wagner et al. - OCD (Classic)
params.standardize = true;
params.combine = true;
[stopGeneration state] = oscTaxonomy('nsga2', 'uf2', 2, 1000, 1, ...
    {'hv', 'epsilon', 'r'}, {'std', 'reg'}, {'adaptTest'}, params)
% Bui et al. - DQP
[stopGeneration state] = oscTaxonomy('nsga2', 'uf2', 2, 1000, 1, ...
    {'dqp'}, {'direct'}, {'threshold'}, params)
% Guerrero et al. - LSSC
params.standardize = false;
params.combine = false;
[stopGeneration state] = oscTaxonomy('nsga2', 'uf2', 2, 1000, 1, ...
    {'mdr', 'hv', 'epsilon'}, {'reg'}, {'validThreshold'}, params)
% Wagner et al. - OCD-HV
[stopGeneration state] = oscTaxonomy('nsga2', 'uf2', 2, 1000, 1, ...
    {'hv'}, {'std'}, {'adaptTest'}, params)
% Goel and Stander - CR
[stopGeneration state] = oscTaxonomy('nsga2', 'uf2', 2, 1000, 1, ...
    {'cr'}, {'moving'}, {'threshold'}, params)