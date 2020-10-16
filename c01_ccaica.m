clear;clf;clc
rootdir = getenv('ROOTDIR');
addpath(fullfile(getenv('FREESURFER_HOME'),'matlab'));
addpath /data/EDB/opt/palm.git;
addpath /data/EDB/opt/toolbox.git/share;
addpath(fullfile(rootdir,'code','robustica'));
warning off all

runnum  = 'run19';
nP       = [1 1]*10000;
C = {'NIMH','CMI'};
P = {'partial'};
netmatsdir = 'netmats';
numMins = 9;
exclude = true;
outprefix = '4to1';
switch outprefix
    case '2to1'
        NpcaY   = [64 134];
        NpcaX   = [29 29];
    case '3to1'
        NpcaY   = [32 79];
        NpcaX   = [29 29];
    case '4to1'
        NpcaY   = [17 53];
        NpcaX   = [29 29];
end

% TRs
TR = cell(size(C));
TR{strcmpi(C,'NIMH')} = [300 2];
TR{strcmpi(C,'CMI')}  = [415 1.45; 740 .8 ; 745 .8; 370 .8];

% Load the structure names
[~,~,Lctab] = read_annotation(fullfile(rootdir,'atlases','Schaefer2018_LocalGlobal','Parcellations','FreeSurfer5.3','fsaverage','label','lh.Schaefer2018_200Parcels_17Networks_order.annot'));
[~,~,Rctab] = read_annotation(fullfile(rootdir,'atlases','Schaefer2018_LocalGlobal','Parcellations','FreeSurfer5.3','fsaverage','label','rh.Schaefer2018_200Parcels_17Networks_order.annot'));
Lnames = Lctab.struct_names(2:end);
Rnames = Rctab.struct_names(2:end);
Snames = {...
    'LH_Accumbens'
    'LH_Amygdala'
    'LH_Caudate'
    'LH_Hippocampus'
    'LH_Pallidum'
    'LH_Putamen'
    'LH_Thalamus'
    'LH_VentralDC'
    'RH_Accumbens'
    'RH_Amygdala'
    'RH_Caudate'
    'RH_Hippocampus'
    'RH_Pallidum'
    'RH_Putamen'
    'RH_Thalamus'
    'RH_VentralDC'};
all_names = vertcat(Lnames,Rnames,Snames);

% Initialise variables for later use:
Y0       = cell(2,1);
X0       = cell(2,1);
Z0       = cell(2,1);
Xnames   = cell(2,1);
Znames   = cell(2,1);
Y1       = cell(2,1);
X1       = cell(2,1);
Z1       = cell(2,1);
Q        = cell(2,1);
Y1       = cell(2,1);  for c = 1:numel(Y1),   Y1{c}     = cell(2,1); end
U        = cell(2,1);  for c = 1:numel(U),     U{c}     = cell(2,1); end
V        = cell(2,1);  for c = 1:numel(V),     V{c}     = cell(2,1); end
R        = cell(2,1);  for c = 1:numel(R),        R{c}        = cell(2,1); end
pvals    = cell(2,1);  for c = 1:numel(pvals),    pvals{c}    = cell(2,1); end
A        = cell(2,1);  for c = 1:numel(A),     A{c}     = cell(2,1); end
B        = cell(2,1);  for c = 1:numel(B),     B{c}     = cell(2,1); end
Aload    = cell(2,1);  for c = 1:numel(Aload),    Aload{c}    = cell(2,1); end
Bload    = cell(2,1);  for c = 1:numel(Bload),    Bload{c}    = cell(2,1); end
Ryx      = cell(2,1);  for c = 1:numel(Ryx),      Ryx{c}      = cell(2,1); end
Rxy      = cell(2,1);  for c = 1:numel(Rxy),      Rxy{c}      = cell(2,1); end
Ryxprop  = cell(2,1);  for c = 1:numel(Ryxprop),  Ryxprop{c}  = cell(2,1); end
Rxyprop  = cell(2,1);  for c = 1:numel(Rxyprop),  Rxyprop{c}  = cell(2,1); end
newnet   = cell(2,1);  for c = 1:numel(newnet),   newnet{c}   = cell(2,1); end
K        = cell(2,1);  for c = 1:numel(K),        K{c}        = cell(2,1); end
Y1avgnet = cell(2,1);  for c = 1:numel(Y1avgnet), Y1avgnet{c} = cell(2,1); end

% For each cohort
for c = 1:numel(C)
    
    % For each residualisation scheme (full/partial)
    for p = 1:numel(P)
        fprintf('%s, %s\n',C{c},P{p});
        
        % Load the variables of interest, and confounds
        X0{c}     = palm_strcsvread(fullfile(rootdir,'stats',runnum,sprintf('%s_X.csv',C{c})));
        Z0{c}     = palm_strcsvread(fullfile(rootdir,'stats',runnum,sprintf('%s_Z.csv',C{c})));
        Xnames{c} = X0{c}(1,2:end);
        Sids      = X0{c}(2:end,1);
        
        % For each subject, load the netmats
        for s = 1:numel(Sids)
            if      exist( fullfile(rootdir,sprintf('derivatives_%s',C{c}),netmatsdir,Sids{s},'func',sprintf('%s_task-rest_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',Sids{s},P{p})),'file')
                % This loads the single, longer run, for the CMI sample
                netfile  = fullfile(rootdir,sprintf('derivatives_%s',C{c}),netmatsdir,Sids{s},'func',sprintf('%s_task-rest_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',Sids{s},P{p}));
                lenfile0 = fullfile(rootdir,sprintf('derivatives_%s',C{c}),'netmats', Sids{s},'func',sprintf('%s_task-rest_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.txt',Sids{s},P{p}));
                lenfile  = fullfile(rootdir,sprintf('derivatives_%s',C{c}),netmatsdir,Sids{s},'func',sprintf('%s_task-rest_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.txt',Sids{s},P{p}));
            elseif  exist( fullfile(rootdir,sprintf('derivatives_%s',C{c}),netmatsdir,Sids{s},'func',sprintf('%s_task-rest_run-both_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',Sids{s},P{p})),'file')
                % This loads the two merged runs, for the CMI sample
                netfile  = fullfile(rootdir,sprintf('derivatives_%s',C{c}),netmatsdir,Sids{s},'func',sprintf('%s_task-rest_run-both_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',Sids{s},P{p}));
                lenfile0 = fullfile(rootdir,sprintf('derivatives_%s',C{c}),'netmats', Sids{s},'func',sprintf('%s_task-rest_run-both_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.txt',Sids{s},P{p}));
                lenfile  = fullfile(rootdir,sprintf('derivatives_%s',C{c}),netmatsdir,Sids{s},'func',sprintf('%s_task-rest_run-both_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.txt',Sids{s},P{p}));
            elseif  exist( fullfile(rootdir,sprintf('derivatives_%s',C{c}),netmatsdir,Sids{s},'ses-1','func',sprintf('%s_ses-1_task-rest_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',Sids{s},P{p})),'file')
                % This is for the NIMH sample (always a single run)
                netfile  = fullfile(rootdir,sprintf('derivatives_%s',C{c}),netmatsdir,Sids{s},'ses-1','func',sprintf('%s_ses-1_task-rest_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',Sids{s},P{p}));
                lenfile0 = fullfile(rootdir,sprintf('derivatives_%s',C{c}),'netmats', Sids{s},'ses-1','func',sprintf('%s_ses-1_task-rest_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.txt',Sids{s},P{p}));
                lenfile  = fullfile(rootdir,sprintf('derivatives_%s',C{c}),netmatsdir,Sids{s},'ses-1','func',sprintf('%s_ses-1_task-rest_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.txt',Sids{s},P{p}));
            else
                warning('File not found for subject: %s', Sids{s});
                netfile = '';
            end
            if ~ isempty(netfile)
                y  = load(netfile);
                n0 = load(lenfile0);
                n  = load(lenfile);
            else
                y  = NaN;
                n0 = NaN;
                n  = NaN;
            end
            if s == 1
                Y0{c} = zeros(numel(Sids),size(y,1));
                nT0   = zeros(numel(Sids),1);
                nT    = zeros(numel(Sids),1);
            end
            Y0{c}(s,:) = y';
            nT0(s,1) = n0;
            nT (s,1) = n;
        end
        X0{c}     = cell2mat(X0{c}(2:end,2:end));
        Znames{c} = Z0{c}(1,2:end);
        Z0{c}     = cell2mat(Z0{c}(2:end,2:end));
        
        % Rescale age and IQ, so as to have nicer ranges, better condition
        % numbers
        idx = find(strcmpi(Znames{c},'Age')|strcmpi(Znames{c},'IQ')|strcmpi(Znames{c},'IQ (WASI)'));
        Z0{c}(:,idx) = Z0{c}(:,idx)/10;
        
        % Remove participants who we wanted to exclude and marked as NaN above.
        idx = isnan(nT);
        X0{c}(idx,:) = [];
        Z0{c}(idx,:) = [];
        Y0{c}(idx,:) = [];
        nT0(idx)     = [];
        nT(idx)      = [];
        
        % Remove participants with incomplete runs
        idx = false(size(nT));
        for r = 1:size(TR{c},1)
            idx = nT0 == TR{c}(r,1) | idx;
        end
        Y0{c}(~idx,:) = [];
        X0{c}(~idx,:) = [];
        Z0{c}(~idx,:) = [];
        nT0  (~idx)   = [];
        nT   (~idx)   = [];
        
        % Remove participants with too few volumes (this only matters if
        % netmats were computed after scrubbing)
        for r = 1:size(TR{c},1)
            idx = nT0 == TR{c}(r,1);
            nT(idx) = nT(idx)*TR{c}(r,2);
        end
        idx = nT < numMins*60;
        Y0{c}(idx,:) = [];
        X0{c}(idx,:) = [];
        Z0{c}(idx,:) = [];
        nT0  (idx)   = [];
        nT   (idx)   = [];
        
        % Drop subjects with missing entries in X or Z
        nE    = size(Y0{c},2); % number of network edges
        idxnotmiss = ~ any(isnan([X0{c} Z0{c}]),2);
        X1{c} = X0{c}(idxnotmiss,:);
        Z1{c} = Z0{c}(idxnotmiss,:);
        Y1{c}{p} = Y0{c}(idxnotmiss,:);

        % Final sample size
        N     = size(X1{c},1);
        fprintf('Sample size: %d\n',N);
        
        % Add an intercept
        Z1{c} = horzcat(Z1{c},ones(size(Z1{c},1),1));

        % Do the CCA
        [pvalstmp,cc,~,~,~,~,Aq,Bq,Uq,Vq,Q{c}] = ...
            permccaica(Y1{c}{p},X1{c},Z1{c},NpcaY(c),NpcaX(c),nP(p));
        pvals{c}{p} = pvalstmp;
        R{c}{p}     = cc;
        A{c}{p}     = Aq;
        B{c}{p}     = Bq;
        U{c}{p}     = Uq;
        V{c}{p}     = Vq;
        
        fprintf('p-values:\n')
        disp(pvals{c}{p})
        K{c}{p} = numel(pvals{c}{p});
        
        % Average (residualized) network, for display purposes (if at all
        % of interest).
        Y1avgnet{c}{p} = rewrap(mean(Q{c}'*Y1{c}{p}));
    end
end
%%
% Scatter plots with the canonical correlations:
for c = 1:numel(C)
    for p = 1:numel(P)
        for k = 1:sum(pvals{c}{p} <= 0.05)
            if sum(pvals{c}{p} <= 0.05) >= 1
            scatter(Q{c}'*U{c}{p}(:,k),Q{c}'*V{c}{p}(:,k),'.')
            print('-dpng','-r100',fullfile(rootdir,'stats',runnum,sprintf('correlation_between_canonical_variables_%s_%s_k%d.png',C{c},P{p},k)));
            end
        end
    end
end

% Within each site, compute the correlation between input and
% canonical variables:
k = 4;
for c = 1:numel(C)
    for p = 1:numel(P)
        subplot(1,2,1);
        Aload{c}{p} = corr(Q{c}'*Y1{c}{p},Q{c}'*U{c}{p});
        imagesc(Aload{c}{p}(:,1:k)); colorbar
        ylabel(sprintf('Original variables %s',C{c}));
        xlabel(sprintf('Canonical variables %s',C{c}));
        title('Imaging side');
        daspect([1 1 1])
        
        subplot(1,2,2)
        Bload{c}{p} = corr(Q{c}'*X1{c},Q{c}'*V{c}{p});
        imagesc(Bload{c}{p}(:,1:k)); colorbar
        title('Clinical side');
        ylabel(sprintf('Original variables %s',C{c}));
        xlabel(sprintf('Canonical variables %s',C{c}));
        daspect([1 1 1])
        if ~palm_isoctave
            sgtitle(sprintf('Loadings, %s, %s',C{c},P{p}))
        end
        print('-dpng','-r100',fullfile(rootdir,'stats',runnum,sprintf('loadings_%s_%s.png',C{c},P{p})));
    end
end

%% 
% If two sites were used, compute the cross-correlations among
% canonical variables
%corr(Q{1}'*Y1{1}*Arec{1}{1},Ures{1}{1}) -> this is identity
%corr(Q{1}'*X1{1}*Brec{1}{1},Vres{1}{1}) -> this is identity
if numel(C) == 2 && (nP(p) >= 100)
    for c = 1:numel(C)
        if c == 1
            d = 2;
        elseif c == 2
            d = 1;
        end
        for p = 1:numel(P)
            nsigc = sum(pvals{c}{p} <= 0.05);
            nsigd = sum(pvals{d}{p} <= 0.05);
            if nsigc >= 1 && nsigd >= 1
                fprintf('Imaging side\n')
                corr(...
                    Q{c}'*Y1{c}{p}*A{c}{p}(:,1:nsigc),...  %Calculate the U for a specific dataset with its own weights
                    Q{c}'*Y1{c}{p}*A{d}{p}(:,1:nsigd))     %Calculate the U using the same dataset using weights from the other dataset
                [pu,pc] = permcorr(...
                    Q{c}'*Y1{c}{p}*A{c}{p}(:,1:nsigc),...
                    Q{c}'*Y1{c}{p}*A{d}{p}(:,1:nsigd),10000)
                [~,~,padj] = fdr(pu(:));
                pfdr = reshape(padj,size(pu))
                [sum(pu(:)<=.05) sum(pc(:)<=.05) sum(pfdr(:)<=.05)]
                fprintf('Clinical side\n')
                corr(...
                    Q{c}'*X1{c}*B{c}{p}(:,1:nsigc),...
                    Q{c}'*X1{c}*B{d}{p}(:,1:nsigd))
                [pu,pc] = permcorr(...
                    Q{c}'*X1{c}*B{c}{p}(:,1:nsigc),...
                    Q{c}'*X1{c}*B{d}{p}(:,1:nsigd),10000)
                [~,~,padj] = fdr(pu(:));
                pfdr = reshape(padj,size(pu))
                [sum(pu(:)<=.05) sum(pc(:)<=.05) sum(pfdr(:)<=.05)]
            end
        end
    end
end

%% Compute the redundancy index
for c = 1:numel(C)
    for p = 1:numel(P)
        Ryx{c}{p}     = R{c}{p}.^2.*mean(Aload{c}{p}.^2,1); % column IV of Stewart & Love (1968)
        Ryxprop{c}{p} = Ryx{c}{p}/sum(Ryx{c}{p});           % column V of Stewart & Love (1968)
        Rxy{c}{p}     = R{c}{p}.^2.*mean(Bload{c}{p}.^2,1); % column IV of Stewart & Love (1968)
        Rxyprop{c}{p} = Rxy{c}{p}/sum(Rxy{c}{p});           % column V of Stewart & Love (1968)
    end
end

%% Wrap back to a netmat
nN = max(roots([1 -1 -2*nE])); % number of network nodes
newnet = cell(2,1); for c = 1:numel(newnet), newnet{c} = cell(2,1); end
for c = 1:numel(C)
    for p = 1:numel(P)
        newnet{c}{p} = zeros(nN,nN,K{c}{p});
        for k = 1:K{c}{p}
            newnet{c}{p}(:,:,k) = rewrap(Aload{c}{p}(:,k));
        end
    end
end

% Save with the nodes in a particularly good order:
node_names = strcsvread(fullfile(rootdir,'code','node_names.csv'));
[~,idx] = sort(cell2mat(node_names(:,2)));
node_names_sorted = node_names(idx,:)';
for c = 1:numel(C)
    mkdir(fullfile(rootdir,'stats','netjs',outprefix,C{c}));
    fid = fopen(fullfile(rootdir,'stats','netjs',outprefix,C{c},'node_names_sorted.csv'),'w');
    fprintf(fid,'%d,%d,%s,%s,%s\n',node_names_sorted{:,:});
    fclose(fid);
end
for c = 1:numel(C)
    for p = 1:numel(P)
        for k = 1:K{c}{p}
            if false%(c == 1 && any(k == [1 3 5])) % Flip some signs
                dlmwrite(fullfile(rootdir,'stats','netjs',outprefix,C{c},sprintf('newnet_c%d_p%d_k%d.csv',c,p,k)),-newnet{c}{p}(idx,idx,k) + eye(size(newnet{c}{p},1)),'delimiter',' ');
            else
                dlmwrite(fullfile(rootdir,'stats','netjs',outprefix,C{c},sprintf('newnet_c%d_p%d_k%d.csv',c,p,k)),+newnet{c}{p}(idx,idx,k) + eye(size(newnet{c}{p},1)),'delimiter',' ');
            end
            if false%(c == 2 && any(k == [1 2 5])) % Flip some signs
                dlmwrite(fullfile(rootdir,'stats','netjs',outprefix,C{c},sprintf('newnet_c%d_p%d_k%d.csv',c,p,k)),-newnet{c}{p}(idx,idx,k) + eye(size(newnet{c}{p},1)),'delimiter',' ');
            else
                dlmwrite(fullfile(rootdir,'stats','netjs',outprefix,C{c},sprintf('newnet_c%d_p%d_k%d.csv',c,p,k)),+newnet{c}{p}(idx,idx,k) + eye(size(newnet{c}{p},1)),'delimiter',' ');
            end
        end
        dlmwrite(fullfile(rootdir,'stats','netjs',outprefix,C{c},sprintf('netmat_c%d_p%d.csv',c,p)),tanh(Y1avgnet{c}{p}(idx,idx)) + eye(size(Y1avgnet{c}{p},1)),'delimiter',' ');
        dlmwrite(fullfile(rootdir,'stats','netjs',outprefix,C{c},sprintf('linkage_c%d_p%d.csv',c,p)),linkage(Y1avgnet{c}{p}(idx,idx) + eye(size(Y1avgnet{c}{p},1))),'delimiter',' ');
    end
end

%%
[~,tmp] = system('whoami');
whoami = deblank(tmp);
save(fullfile(rootdir,'stats',runnum,sprintf('%s_%s_%d_%s_%s_%s_CCA+ICA_partial+full.mat',outprefix,whoami,max(nP),version('-release'),lower(computer),date)))