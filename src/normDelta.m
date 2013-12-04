function [ resA, resB] = normDelta( setA, setB )
%NORMDELTA Returns the number or elements in the sets determined by the
%delta function.
%
%   The delta function returns the set of the elements of setA that are
%   dominated by at least one element of setB and vice versa. 
%
%   This function returns the relative number of elements in those sets.
%
%   It asumes that the elements of setA and setB are non-dominated.
%
%   Parameters:
%   * setA, a matrix of nA rows and m columns.
%   * setB, a matrix of nB rows and m columns.
%
%   Results:
%   * resA, ratio of elements of setA in delta(setA, setB).
%   * resA, ratio of elements of setB in delta(setB, setA).
%
%   See:
%   Marti, L., Garcia, J., Berlanga, A., & Molina, J. M. (2009).
%   An Approach to Stopping Criteria for Multi-Objective Optimization
%   Evolutionary Algorithms: The MGBM Criterion. In 2009 IEEE Conference
%   on Evolutionary Computation (CEC 2009), Piscataway, New Jersey
%
%   * Could be improved by using:
%   http://www.mathworks.cn/matlabcentral/fileexchange/17251
%
% ---
% Version: $Id: normDelta.m 97 2010-10-09 16:33:46Z wagner $

nA = size(setA,1);
nB = size(setB,1);

parIndex = paretofront([setA; setB]);
resA = sum(~parIndex(1:nA));
resB = sum(~parIndex(nA+1:end));

resA = resA/nA;
resB = resB/nB;