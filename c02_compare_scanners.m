
%rootdir = getenv('/gpfs/gsfs6/users/EDB/MErest');
%addpath('/gpfs/gsfs6/users/EDB/opt/toolbox.git/share');
%addpath('/gpfs/gsfs6/users/EDB/opt/palm.git');


%load(fullfile(rootdir,sprintf('stats/run19/2to1_winkleram_10000_2020a_glnxa64_18-Jun-2020_CCA+ICA_partial+full.mat')))

%load(fullfile(rootdir,sprintf('stats/run19/2to1_winkleram_10000_2020a_glnxa64_03-Jul-2020_CCA+ICA_partial+full.mat')))


% For each cohort
KSstat = cell(size(C));
KSpval = cell(size(C));
KSpfwe = cell(size(C));
KSprint = cell(size(C));
for c = 1:numel(C)
    KSstat{c} = cell(size(P));
    KSpval{c} = cell(size(P));
    KSpfwe{c} = cell(size(P));
    KSprint{c} = cell(size(P));
end

for c = 1:numel(C)
    for p = 1:numel(P)
        fprintf('Working on %s %s\n',C{c},P{p});
        idxa = contains(Znames{c},'Scanner');
        idxb = contains(Znames{c},'Site');
        idx  = find(idxa | idxb);
        [~,~,ScannerIdx] = unique(Z1{c}(:,idx),'rows');
        nScanners = numel(unique(ScannerIdx));
        nTests = nScanners*(nScanners-1)/2;
        KSstat{c}{p} = zeros(K{c}{p},nTests);
        KSpval{c}{p} = zeros(K{c}{p},nTests);
        KSpfwe{c}{p} = zeros(K{c}{p},nTests);
        sc = 1;
        for scanner1 = 1:nScanners-1
            for scanner2 = scanner1+1:nScanners
                fprintf('- Testing Scanner %d vs Scanner %d\n',scanner1,scanner2);
                Uscanner1 = U{c}{p}(ScannerIdx == scanner1,:);
                Uscanner2 = U{c}{p}(ScannerIdx == scanner2,:);
                for k = 1:size(Uscanner1,2)
                    [~,KSpval{c}{p}(k,sc),KSstat{c}{p}(k,sc)] = kstest2(Uscanner1(:,k),Uscanner2(:,k));
                end
                sc = sc + 1;
            end
        end
        KSpfwe{c}{p} = 1-(1-KSpval{c}{p}).^(nTests*K{c}{p});
    end
end

for c = 1:numel(C)
    for p = 1
        tmp = cat(3,KSstat{c}{p},KSpval{c}{p},KSpfwe{c}{p});
        tmp = permute(tmp,[3 1 2]);
        tmp = reshape(tmp,K{c}{p}*3,size(KSstat{c}{p},2));
        KSprint{c}{p} = tmp;
    end
end
