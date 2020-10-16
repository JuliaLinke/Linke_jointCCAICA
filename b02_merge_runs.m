rootdir = getenv('ROOTDIR');
C = {'full', 'partial'};
tmp = dir(fullfile(rootdir,'derivatives_CMI','netmats','sub-*'));
S = cell(numel(tmp),1);
for s = 1:numel(tmp)
    S{s} = tmp(s).name;
end
for c = 1:numel(C)
    for s = 1:numel(S)
        if      exist(fullfile(rootdir,'derivatives_CMI','netmats',S{s},'func',sprintf('%s_task-rest_run-1_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',S{s},C{c})),'file') && ...
                exist(fullfile(rootdir,'derivatives_CMI','netmats',S{s},'func',sprintf('%s_task-rest_run-2_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',S{s},C{c})),'file')
            fprintf('%s, %s, two runs\n',S{s},C{c})
            r1 = load(fullfile(rootdir,'derivatives_CMI','netmats',S{s},'func',sprintf('%s_task-rest_run-1_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',S{s},C{c})));
            r2 = load(fullfile(rootdir,'derivatives_CMI','netmats',S{s},'func',sprintf('%s_task-rest_run-2_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',S{s},C{c})));
            n1 = load(fullfile(rootdir,'derivatives_CMI','netmats',S{s},'func',sprintf('%s_task-rest_run-1_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.txt',S{s},C{c})));
            n2 = load(fullfile(rootdir,'derivatives_CMI','netmats',S{s},'func',sprintf('%s_task-rest_run-2_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.txt',S{s},C{c})));
            r = (r1 + r2)/2*sqrt(2);
            n = n1 + n2;
            dlmwrite(fullfile(rootdir,'derivatives_CMI','netmats',S{s},'func',sprintf('%s_task-rest_run-both_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',S{s},C{c})),r,'precision','%0.18e');
            dlmwrite(fullfile(rootdir,'derivatives_CMI','netmats',S{s},'func',sprintf('%s_task-rest_run-both_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.txt',S{s},C{c})),n,'precision','%d');
        elseif  exist(fullfile(rootdir,'derivatives_CMI','netmats',S{s},'func',sprintf('%s_task-rest_netmat-%s_atlas-Schaefer2018-200P+17N_space-T1w.csv',S{s},C{c})),'file')
            fprintf('%s, %s, one run\n',S{s},C{c})
        else
            fprintf('%s, %s, no usable runs\n',S{s},C{c});
        end
    end
end
