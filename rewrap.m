function W = rewrap(V)
nE = length(V);
nN  = max(roots([1 -1 -2*nE])); % number of network nodes
idx = logical(triu(ones(nN),+1));
tmp = zeros(nN,nN);
tmp(idx) = V;
W = tmp + tmp';


