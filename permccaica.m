function [pvals,cc1,Aqp,Bqp,Uqp,Vqp,Aq,Bq,Uq,Vq,Q] = permccaica(Y,X,Z,NpcaY,NpcaX,nP)
% X, Y, Z        = input variables
% NpcaX, NpcaY   = Number of PCA components to retain
% nP             = Number of  permutations

% Outputs:
% pvals          = FWER-corrected pvals using closure
% cc1            = Canonical correlations post ICA (unpermuted)
% Aqp, Bqp, Uqp, Vqp = Canonical coefficients and variables in the PCA/Q
%                  space (i.e., N' by NpcaY and N' by NpcaX).
% Aq, Bq, Uq, Vq = Canonical coefficients and variables in the original
%                  spaces (i.e., N by Ny and N by Nx).
% Q              = Semi-orthogonal residualisation matrix.s


% Deconfound using Huh-Jhun:
N     = size(X,1);
[Q,D] = svd(null(Z'));
Q     = Q*D;

% Remove nuisance, deconfound as defined above:
Yq    = Q'*Y;      % Y residualised and exchangeable, N' by Ny
Xq    = Q'*X;      % X residualised and exchangeable, N' by Nx
Nnew  = size(Q,2); % same as N' = N - rank(Z) if Z is full rank
df    = Nnew - 1;  % used to scale CCA outputs to have unit variance

% PCA for dim reduction if requested
[Yqp,LY,~,extras] = epca(Yq,rank(Yq)); % Yq = Yqp*LY'
explainedY = extras.explained;
nCompY = min(NpcaY,size(Yqp,2));
if NpcaY > 1
    Yqp    = Yqp(:,1:nCompY);
    LY     = LY (:,1:nCompY);
elseif NpcaY < 1  % if the number of requested PCs is <1, use it as fraction of variance explained
    idx    = cumsum(explainedY') <= NpcaY;
    nCompY = sum(idx);
    Yqp    = Yqp(:,idx);
    LY     = LY(:,idx);
end
[Xqp,LX,~,extras] = epca(Xq,rank(Xq)); % Xq = Xqp*LX'
explainedX = extras.explained;
nCompX = min(NpcaX,size(Xqp,2));
if NpcaX > 1
    Xqp    = Xqp(:,1:nCompX);
    LX     = LX (:,1:nCompX);
elseif NpcaX < 1 % if the number of requested PCs is <1, use it as fraction of variance explained
    idx    = cumsum(explainedX') <= NpcaX;
    nCompX = sum(idx);
    Xqp    = Xqp(:,idx);
    LX     = LX(:,idx);
end
fprintf('Number of variables [img cln]:\n')
disp([nCompY nCompX]);
fprintf('Variance explained on the imaging side (%%):\n');
disp(cumsum(explainedY'*100));
fprintf('Variance explained on the clinical side (%%):\n');
disp(cumsum(explainedX'*100));

% For each permutation
tic
for p = 1:nP
    fprintf('Permutation %d/%d: ',p,nP);
    if p == 1
        % First permutation is no permutation
        idxY = (1:Nnew);
        idxX = (1:Nnew);
    else
        % If not first permutation, permute randomly
        idxX = randperm(Nnew);
        idxY = randperm(Nnew);
    end
    
    % This is the "sequential" way
    % For each canonical component, do the CCA from k onwards
    % This will not affect the first permutation, but will effectively
    % remove from later k the earlier ones. Since these later ones will
    % only be tested if the earlier ones were significant (step-down)
    % this means the earlier ones will act as nuisance, and will not be
    % permuted together.
    [ccs,Aqps,Bqps,~,~] = intraperm(Q*Yqp(idxY,:),Q*Xqp(idxX,:),df,'stack');
    Ts = wilks(ccs); % test statistic
    
    % Keep the results of the 1st permutation
    if p == 1 
        cc1   = ccs;
        T1    = Ts;
        Aqp   = Aqps;
        Bqp   = Bqps;
        pvals = zeros(size(ccs));
    end
    
    % Increment counter for the p-value
    pvals = pvals + (T1 <= Ts);
end
fprintf('Time per permutation: %g\n', toc/nP);

% Compute the p-values
pvals = cummax(pvals/nP,2);

% Canonical variables from residualised data after PCA:
Uqp = Q*Yqp*Aqp;
Vqp = Q*Xqp*Bqp;

% Canonical coefficients/variables that refer to residualised data before PCA:
Aq = pinv(LY')*Aqp;
Bq = pinv(LX')*Bqp;
Uq  = Q*Yq*Aq;
Vq  = Q*Xq*Bq;
end

% =================================================================
function [cc,A,B,U,V] = cca(Y,X,df)
% The usual CCA.
[Qy,Ry,Ty] = qr(Y,0);
[Qx,Rx,Tx] = qr(X,0);
K  = min(rank(Y),rank(X));
[L,D,M] = svds(Qy'*Qx,K);
cc = min(max(diag(D(:,1:K))',0),1);
A  = Ry\L(:,1:K)*sqrt(df);
B  = Rx\M(:,1:K)*sqrt(df);
A(Ty,:) = A;
B(Tx,:) = B;
U       = Y*A;
V       = X*B;
end

% =================================================================
function [cc,A,B,U,V] = ccaica(Y,X,df,method)
% Given two sets of variables, run an initial CCA, producing left 
% and right canonical variables; these are combined using one of two
% methods: 'stack' or 'sum'. Then do an ICA, and combine the canonical
% loadings with the unmixing coefficients. Return unmixing coefficients
% (A and B) that act on the initial inputs, the separated
% sources (U and V), and the diagonal of their correlation (cc).
% Sources and coefficients are sorted such thay corr(U,V) has a
% dominant diagonal.
[~,Amix,Bmix,Umix,Vmix] = cca(Y,X,df);
switch method
    case 'sum'
        UVmix = Umix+Vmix;
    case 'stack'
        UVmix = vertcat(Umix,Vmix);
end
[~,mixUV,~,~] = robustica(UVmix',{'prewhi',false});
unmixUV = inv(mixUV');
A = Amix*unmixUV*sqrt(df/(df+1)); %#ok<MINV>
B = Bmix*unmixUV*sqrt(df/(df+1)); %#ok<MINV>
U = Y*A;
V = X*B;
[cc,idxU,idxV,signflip] = todiag(bcorr(U,V));
cc = diag(cc)';
A = A(:,idxU)*signflip;
B = B(:,idxV);
U = U(:,idxU)*signflip;
V = V(:,idxV);
end

% =================================================================
function [cc,A,B,U,V] = icapair(Umix,Vmix,df,method) % NOT IN USE.
% Given two sets of prewhitened variables (such as after CCA),
% combine them using one of two methods: 'stack' or 'sum'.
% do an ICA, and return unmixing coefficients (A and B), separated
% sources (U and V), and the diagonal of their correlation (cc).
% Sources and coefficients are sorted such thay corr(U,V) has a
% dominant diagonal.
% The main difference between this and "ccaica" is that this one
% assumes inputs are not subjected to CCA again everytime (i.e. for
% every k). Probably a bad idea, and not currently in use.
switch method
    case 'sum'
        UVmix = Umix+Vmix;
    case 'stack'
        UVmix = vertcat(Umix,Vmix);
end
[~,mixUV,~,~] = robustica(UVmix',{'prewhi',false});
unmixUV = inv(mixUV')*sqrt(df/(df+1)); %#ok<MINV>
U  = Y*unmixUV;
V  = X*unmixUV;
[cc,idxU,idxV,signflip] = todiag(bcorr(U,V));
cc = diag(cc)';
A  = unmixUV(:,idxU)*signflip;
B  = unmixUV(:,idxV);
U  = U(:,idxU)*signflip;
V  = V(:,idxV);
end

% =================================================================
function [cc,A,B,U,V] = intraperm(Y,X,df,method)
% This is the main function, that calls ccaica, cumulatively do the
% products for the coefficients such that they can be applied to the
% orginal data to produce canonical post-ICA variables, diagonalises
% via permutation of rows and columns the correlation matrix, sort
% them in descending order, and flips signs to make sure all
% correlations are positive.
Yk = Y;
Xk = X;
K  = min(rank(Y),rank(X));
ccs = zeros(size(K,1));
for k = 1:K
    fprintf('%d ',k)
    [cck,Ak,Bk] = ccaica(Yk,Xk,df,method);
    ccs(k) = cck(1);
    A = [Ak null(Ak')];
    B = [Bk null(Bk')];
    if k == 1
        Aks = A;
        Bks = B;
    else
        Aks(:,k:end) = Aks(:,k:end)*A;
        Bks(:,k:end) = Bks(:,k:end)*B;
    end
    Yk = Yk*A(:,2:end);
    Xk = Xk*B(:,2:end);
end
fprintf('\n');
Uks = Y*Aks;
Vks = X*Bks;
[cc,idxU,idxV,signflip] = todiag(bcorr(Uks,Vks));
cc = diag(cc)';
A = Aks(:,idxU(1:K))*signflip;
B = Bks(:,idxV(1:K));
U = Uks(:,idxU(1:K))*signflip;
V = Vks(:,idxV(1:K));
end

% =================================================================
function [Y,idxR,idxC,S] = todiag(X)
% Takes a matrix X and ensures it has a dominant diagonal based
% on largest values. The reordering of rows and columns is determined
% by the returned indices idxR and idxC: Y = X(idxR,idxC). Return also
% a sign-flipping matrix that can be used to swap signs of the original
% variables used to create X if X is a correlation or covariance matrix.
Sx      = sign(X); % sign is preserved for later
X       = abs(X);  % ordering is based on absolute values
Y       = X;
[nR,nC] = size(X);
N       = min(nR,nC);
idxR    = zeros(1,N);
idxC    = zeros(1,N);
mi      = min(Y(:));
for n = 1:N
    ma      = max(Y(:));
    [r,c]   = find(Y == ma);
    idxR(n) = r;
    idxC(n) = c;
    Y(r,:)  = mi;
    Y(:,c)  = mi;
end
idxR  = [idxR setdiff(1:nR,idxR)];
idxC  = [idxC setdiff(1:nC,idxC)];
Y     = X(idxR,idxC).*Sx(idxR,idxC);
S     = diag(sign(diag(Y)));
end

% =================================================================
function C = bcorr(X,Y)
% Basic equivalent to "corr" that works in Octave and in Matlab
% versions that don't have broadcasting.
X = bsxfun(@minus,  X,mean(X,1));
Y = bsxfun(@minus,  Y,mean(Y,1));
X = bsxfun(@rdivide,X,std(X,1));
Y = bsxfun(@rdivide,Y,std(Y,1));
C = (X'*Y)/size(X,1);
end

% =================================================================
function lW = wilks(cc)
% Compute Wilks' lambda test statistic.
lW = -fliplr(cumsum(fliplr(log(1-cc.^2))));
end

% =================================================================
function R = roy(cc) % NOT IN USE
% Compute Roy's test statistic.
R = cc.^2;
end

% =================================================================
function T = pillai(cc) % NOT IN USE
% Compute Pillai's trace statistic.
T = fliplr(cumsum(fliplr(cc)));
end