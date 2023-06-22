% 0) data preparation: 00_ming_rds_to_txt.R
% mv datasets_list.txt and geneid.txt to 'sn_brain_100/'


% 1) run C3D
% clc
% clear
% cd sn_brain_100/
% file_path=pwd
% cd ..
% ming_C3D_V.m: ming_C3D_V(file_path, 0.01, 1000, 1, 'abs') %35min for 60nets, 40min for 100nets

% 2) run nodeSelection_fdrtool.R

% 3) 
if(0)
    % if use Z-score to select node for each module
    x=load('sc_brain_58/vx_lambda.mat')
    %x=load('sn_brain_100/vx_lambda.mat')
    V=x.V;
    xita=2.5
    [ candidateModules ] = moduleNodesSelection( V, xita )
end

% if use fdrtool qvalue to select node
% read in modules info
%tmp=load('sc_candidateModules.mat')
%tmp=load('sn_candidateModules.mat')
%tmp=load('sc_female_candi.modules.mat')
%tmp=load('sc_female_candi.modules.mat')
%tmp=load('sn_female_candi.modules.mat')
%tmp=load('sn_male_candi.modules.mat')
%tmp=load('sn_male_candi.modules_rep2.mat')
tmp=load('sn_500cells_modules.mat')
%tmp=load('sn_500cells_modules.mat')
%tmp=load('sc_500cells_modules.mat')
modules=tmp.module_df;

N=max(modules.id);
candidateModules = cell(N, 1);
for i = 1:N
  candidateModules{i}=modules.gene(modules.id==i)
end
modulesFinal = candidateModules;

% begin module merging
if(0)
    HPI = setSimilarity( candidateModules );
    size(HPI)
    modulesFinal = candidateModules;

    for i = 1:size(HPI, 1)-1
        for j = (i+1):size(HPI, 2)
            if HPI(i,j)>0.5 % merge these two modules
                [Y, I] = max([length(candidateModules{i}), length(candidateModules{j})]); %big.value, big.index
                if I == 1
                    modulesFinal{j} = [];                
                    HPI(j, :) = zeros(1, size(HPI, 2));
                    HPI(:, j) = zeros(size(HPI, 1), 1);
                else
                    modulesFinal{i} = [];               
                    HPI(i, :) = zeros(1, size(HPI, 2));
                    HPI(:, i) = zeros(size(HPI, 1), 1);
                end
            end
        end
    end

    modulesFinal
end


% Only modules with no less than 5 nodes are kept
i = 1;
while i ~= length(modulesFinal)+1
    if isempty(modulesFinal{i}) || (length(modulesFinal{i})<5)
        modulesFinal(i) = [];
        i = i - 1;
    end
    i = i + 1;
end

modulesFinal

% 4)  module detection across datasets
%% read in multiNetworks
%dinfo = dir('../brain_scRNA-seq/pearson_coexpr.mat_common/*.mat');
dinfo = dir('../brain_snRNA-seq/pearson_coexpr.mat_common/*.mat');
%dinfo = dir('../brain_scRNA-seq/pearson_coexpr.mat_common/female*.mat');
%dinfo = dir('../brain_scRNA-seq/pearson_coexpr.mat_common/male*.mat');
%dinfo = dir('../brain_snRNA-seq/pearson_coexpr.mat_common/female*.mat');
%dinfo = dir('../brain_snRNA-seq/pearson_coexpr.mat_common/male*.mat');
length(dinfo) %41 files
dataset = cell(length(dinfo), 1);
for K = 1 : length(dinfo)
  filename = dinfo(K).name;  %just the name  
  path=dinfo(K).folder;
  dataset{K} = load([path,'/',filename]);
end
length(dataset)

for K = 1:length(dataset)
  %multiNetworks{K}=cell2mat(struct2cell(dataset{K}));
  multiNetworks{K}=abs(cell2mat(struct2cell(dataset{K})));
end
size(multiNetworks{1})


%% V: Module validation
disp('Validation...') %significantModules( modules, multiNetworks, N, permu_times )
num_Nodes= size(multiNetworks{1},1);
tic
[ pvalues_modulePerNet, FDR2 ] = significantModules(modulesFinal, multiNetworks, num_Nodes,100);
runtime = toc;
disp(['Running time: ', num2str(runtime), ' sec.'])


%% save files: Hc, modules_final,pvalues_modulePerNet,FDR2
filenames={dinfo.name};
%save('sc_C3D_out.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
%save('sn_C3D_out.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
%save('sc_C3D_female_out.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
%save('sc_C3D_male_out.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
%save('sn_C3D_female_out.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
%save('sn_C3D_male_out.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
%save('sn_C3D_male_out_rep2.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
%save('sn_C3D_1k.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
save('sn_500cells_modules.fdr.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
%save('sn_1kcells_modules.fdr.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')
%save('sc_500cells_modules.fdr.mat','filenames','modulesFinal','pvalues_modulePerNet','FDR2')

% 5)  run R, 05_GOenrich_plot.R


