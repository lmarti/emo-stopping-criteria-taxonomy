% ========================================================================
% This program calculates Hansen and Jaszkiewicz's R1_R, R2_R and R3_R
% measures for evaluating the quality of sets of nondominated vectors with
% respect to a reference set
%
% Usage:
%   rIndicator(dataFile, referenceSet, idealPoint, nadirPoint, rho, ...
%              method, variant)
%     (all objectives are to be minimized)
%
%   dataFile:     Specifies a matrix (n x d) that contains the output of
%                 a run of a multiobjective optimizer (required)
%
%   referenceSet: Is the name of a matrix (m x d) that contains the
%                 reference set (required)
%
%   idealPoint:   Vector of length d giving the co-ordinates of an ideal
%                 reference point in the objective space (required)
%
%   nadirPoints:  Vector of length d giving the co-ordinates of a nadir
%                 (or worst) point in the objective space. When comparing
%                 two sets of results, it is VERY IMPORTANT that the same
%                 reference point is used. The reference point must be
%                 worse in all dimensions than the worst value of any point
%                 in each dimension (required)
%
%   rho:          This scalar is used in the calculation of the augmented
%                 Tchebycheff function. It controls the degree to which the
%                 Thebycheff function is augmented by a weighted sum
%                 function. Values in the range 0.001 to 0.1 are normal.
%                 This number is ignored if the augmented Tchebycheff
%                 function is not being used (optional, default 0.01)
%
%   method:       Integer in {1,2,3} which specifies whether the
%                 weighted sum (1), Thebycheff (min-max) (2), or
%                 augmented Thebycheff function (3) is used in the
%                 calculation of the R1 and R2 metrics
%                 (optional, default 3)
%
%  variant:       Integer in {1,2,3} which specifieswhich of the R_R
%                 metrics is calculated (optional, default 2)
%
% FURTHER INFORMATION:
% ====================
%    @techreport{hansen98evaluating,
%      author = "Michael Pilegaard Hansen and Andrzej Jaszkiewicz",
%      title = "Evaluating the quality of approximations to the non-dominated set",
%      number = "IMM-REP-1998-7",
%      year = "1998",
%      url = "citeseer.nj.nec.com/hansen98evaluating.html"}
%
%    These measures were reviewed in the PhD thesis of Joshua Knowles and also in
%    a paper by Knowles and Corne published in the 2002 Congress on
%    Evolutionary Computation. These can
%    both be downloaded from http://dbk.ch.umist.ac.uk/knowles/
%    or from the EMOO webpage hosted by Carlos Coello Coello.
% ========================================================================
%
% ---
% Version: $Id: rIndicator.m 33 2010-09-10 15:56:31Z marti-big $

function value = rIndicator(dataFile, referenceSet, idealPoint, ...
    nadirPoint, rho, method, variant)
% check input
if nargin < 4
    error('dataFile, referenceSet, idealPoint, and nadirPoint are required');
elseif nargin < 5
    rho = 0.01;
    method = 3;
    variant = 2;
elseif nargin < 6
    method = 3;
    variant = 2;
elseif nargin < 7
    variant = 2;
    if (method ~= 1 ) && (method ~= 2 ) && (method ~= 3)
        error('method must be 1, 2, or 3');
    end;
end;
if ((method ~= 1 ) && (method ~= 2 ) && (method ~= 3)) || ...
        ((variant ~= 1 ) && (variant ~= 2 ) && (variant ~= 3))
    error('method and variant must be 1, 2, or 3');
end;
if sum(min([dataFile; referenceSet]) < idealPoint) || ...
        sum(max([dataFile; referenceSet]) > nadirPoint)
    error('idealPoint must be smaller and nadirPoint must be greater than all points in dataFile and referenceSet');
end;
if (rho <= 0) || (rho > 1)
    error('rho must be in ]0, 1]');
end;
[n d] = size(dataFile);
[m dm] = size(referenceSet);
if (d < 1 ) || (d ~= dm)
    error('dimensions of dataFile or referenceSet must be equal and > 1');
end;
switch d
    case 2
        weights = load('weights2D.mat');
        weights = weights.ans;
    case 3
        weights = load('weights3D.mat');
        weights = weights.ans;
    case 4
        weights = load('weights4D.mat');
        weights = weights.ans;
    case 5
        weights = load('weights5D.mat');
        weights = weights.ans;
    otherwise        
        weights = createVectors(3, d);
end;
N = size(weights, 1);
maxUtilityData = zeros(N,1);
maxUtilityRef = zeros(N,1);
for i = 1:N
    switch method
        case 1
            utilityData = weightedSum(weights(i,:), dataFile, ...
                idealPoint, nadirPoint);
            utilityRef = weightedSum(weights(i,:), referenceSet, ...
                idealPoint, nadirPoint);
        case 2
            utilityData = regularTchebycheff(weights(i,:), dataFile, ...
                idealPoint, nadirPoint);
            utilityRef = regularTchebycheff(weights(i,:), referenceSet, ...
                idealPoint, nadirPoint);
        case 3
            utilityData = augmentedTchebycheff(weights(i,:), dataFile, ...
                idealPoint, nadirPoint, rho);
            utilityRef = augmentedTchebycheff(weights(i,:), referenceSet, ...
                idealPoint, nadirPoint, rho);            
    end;
    maxUtilityData(i) = max(utilityData);
    maxUtilityRef(i) = max(utilityRef);
end;
switch variant
    case 1
        % R1_R gives the fraction of utility functions on which the
        % reference set is better than the data set. 
        % 0 is the best result. 1 is the worst. 
        % Note: if all approximation sets score 0 then this does not imply 
        % that they are equally good; the reference set is too poor. 
        % Analogously, if all approximations score 1, the reference set is 
        % too good. Try using R2_R or R3_R.
        value = mean( (maxUtilityData < maxUtilityRef) + ...
            0.5.*(maxUtilityData == maxUtilityRef) );        
    case 2
        % R2_R gives the average of the utility difference between the 
        % reference set and the data file. Lower values are better, 
        % with -1 being best and +1 being worst.
        value = mean(maxUtilityRef - maxUtilityData);
    case 3
        % R3_R gives the average of the utility ratio of the data file
        % compared to the reference set. Lower values are better, 
        % with -inf being best and +inf being worst.
        value = mean( (maxUtilityRef - maxUtilityData) ./ ...
            (maxUtilityRef + eps) );
end;

function weights = createVectors(s, d)
c = 1; 
num = nchoosek(s+d-1,d-1);
weights = zeros(num, d);
i = 1;
while i <= (s+1).^d    
    count = int2kary(i,s+1,d);
    loopSum = sum(count);
	if loopSum == s	
        weights(c,:) = count./s;        
        c = c+1;
    end;
    i = i+1;
end;

function kary = int2kary(x, basek, digits)
val = digits-1;
kary = zeros(1, digits);
i = 1;
while x
    if x >= basek.^val
        kary(i) = kary(i)+1;
	    x = x - basek.^val;
    else
        val = val-1;
	    i = i+1;
    end;
end;

function value = weightedSum(weights, objectives, ideal, nadir)
d = length(nadir);
value = ones(size(objectives,1),1);
for i = 1:d
    value = value - weights(i).*( (objectives(:,i) - ideal(i)) ./ ...
        (nadir(i) - ideal(i)) );
end;

function value = regularTchebycheff(weights, objectives, ideal, nadir)
d = length(nadir);
all = zeros(size(objectives));
for i = 1:d
    all(:,d) = weights(i).*( (objectives(:,i) - ideal(i)) ./ ...
        (nadir(i) - ideal(i)) );
end;
value = ones(size(objectives, 1), 1) - max(all, [], 2);

function value = augmentedTchebycheff(weights, objectives, ideal, ...
    nadir, rho)
value = regularTchebycheff(weights, objectives, ideal, nadir) ...
    + rho.*weightedSum(weights, objectives, ideal, nadir);
