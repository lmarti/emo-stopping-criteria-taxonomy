function value = epsilonIndicator(dataFile, referenceSet, method)
%  =======================================================================
%  Implements the unary epsilon indicator as proposed in
%  Zitzler, E., Thiele, L., Laumanns, M., Fonseca, C., and Grunert da
%  Fonseca, V (2003): Performance Assessment of Multiobjective Optimizers:
%  An Analysis and Review. IEEE Transactions on Evolutionary Computation,
%  7(2), 117-132.
%
%  Usage:
%    epsilonIndicator(dataFile, referenceSet, [method])
%      (all objectives are to be minimized)
%
%    datafile:     Specifies a matrix (n x d) that contains the output of
%                  a run of a multiobjective optimizer (required)
%
%    referenceSet: Is the name of a matrix (m x d) that contains the
%                  reference set according to which the indicator is
%                  calculated (required)
%
%    method <0|1>: Determines whether the additive (0) or the
%                  multiplicative (1) version of the indicator is used.
%                  (optional, default 0 (additive epsilon indicator))
%
%  Author:
%    Eckart Zitzler, February 3, 2005 / last update August 9, 2005
%    Tobias Wagner, September 22, 2008 / MATLAB Conversion
%
% ---
% Version: $Id: epsilonIndicator.m 33 2010-09-10 15:56:31Z marti-big $


% check input
if nargin < 2
    error('dataFile and referenceSet are required');
elseif nargin < 3
    method = 0; % unary epsilon indicator
end;
[n d] = size(dataFile);
[m dm] = size(referenceSet);
if (d < 1 ) || (d ~= dm)
    error('dimensions of dataFile or referenceSet must be equal and > 1');
end;
if method == 0
    value = -inf;
else
    value = 0;
end;
for i = 1:m
    for j = 1:n
        for k = 1:d
            switch method
                case 0
                    eps_temp = dataFile(j,k) - referenceSet(i,k);
                case 1
                    if ((referenceSet(i,k)<0)&&(dataFile(j,k)>0)) || ...
                            ((referenceSet(i,k)>0)&&(dataFile(j,k)<0)) ...
                            || (referenceSet(i,k)==0) || (dataFile(j,k)==0)
                        error('dataFile and referenceSet have to be > 0');
                    end;
                    eps_temp = dataFile(j,k) ./ referenceSet(i,k);
            end;
            if k == 1;
                eps_k = eps_temp;
            elseif eps_k < eps_temp
                eps_k = eps_temp;
            end;
        end;
        if j == 1
            eps_j = eps_k;
        elseif eps_j > eps_k
            eps_j = eps_k;
        end;
    end;
    if i == 1
        value = eps_j;
    elseif value < eps_j
        value = eps_j;
    end;
end;
