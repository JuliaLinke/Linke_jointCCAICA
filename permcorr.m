function [P,Pfwer] = permcorr(A,B,nP)
% A and B are two generic matrices with same number of rows, 
% and possibly different number of columns.
% The test will investigate the correlations of all cols of one
% matrix vs all cols of the other matrix.
% nP is the number of permutations
% Output P are the uncorrected p-vals
% Output Pfwer are the corrected p-vals.

[n,ma] = size(A);
[~,mb] = size(B);

for p = 1:nP
    if p == 1
        idxa1 = (1:n)';
        idxa2 = (1:ma)';
        idxb2 = (1:mb)';
    else
        idxa1 = randperm(n)';
        idxa2 = randperm(ma)';
        idxb2 = randperm(mb)';
    end
    T = abs(corr(A(idxa1,idxa2),B(:,idxb2)));
    if p == 1
        C   = zeros(size(T));
        Cmx = C;
        T0  = T;
    end
    C   = C   + (T>=T0);
    Cmx = Cmx + (max(T(:))>=T0);
end
P     = C/nP;
Pfwer = Cmx/nP;