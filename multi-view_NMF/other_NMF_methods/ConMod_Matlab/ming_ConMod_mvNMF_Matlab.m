%% change folder to ConMod_Matlab
%% list all mat files

% S=load('../../propr_coexpr.mat/male_TmY14_perb.mat')
dinfo = dir('../brain_scRNA-seq/pearson_coexpr.mat_common/*.mat');
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



% Calculting the feature networks
tic
disp('Calculating the strengh matrix and the uniformity matrix...')
num_Nodes = size(multiNetworks{1},1);
[Strength, Participation] = featureNets(multiNetworks, num_Nodes);
runtime = toc;
disp(['Running time: ', num2str(runtime), ' sec.'])

%Strength=dlmread('../../strength.txt');
%Participation=dlmread('../../participation.txt');
size(Strength)
size(Participation)

% load information entroyp based integrated network
x=load('../brain_scRNA-seq/integrate.mat')
Information=x.integrate_net;
size(Information)

% Obtaining the candidate modules by multi-view NMF
disp('Obtaining candidate modules by multi-view NMF...')

K = 50;
disp(['K=', num2str(K)])
%X = {Strength, Participation};
%X = {Information, Participation};
X = {Strength, Information};
lambda = [0.01, 0.05];

% multiViewNMF( X, K, lambda, maxIter )
% ~1hr for 6k by 6k matrix, for 1k by 1k matrix
tic
[ H, Hc, objValue ] = multiViewNMF( X, K, lambda, 50 ); 
runtime = toc;
disp(['Running time: ', num2str(runtime), ' sec.'])


% Selecting nodes from the consensus factors
tic
xita = 2; %this output some conserved modules, xita=3 didn't work
modules_final = moduleNodesSelection( Hc, xita );
runtime = toc;
disp(['Running time: ', num2str(runtime), ' sec.'])


%% V: Module validation
disp('Validation...') %significantModules( modules, multiNetworks, N, permu_times )
tic
[ pvalues_modulePerNet, FDR2 ] = significantModules(modules_final, multiNetworks, num_Nodes,100);
runtime = toc;
disp(['Running time: ', num2str(runtime), ' sec.'])


%% save files: Hc, modules_final,pvalues_modulePerNet,FDR2
filenames={dinfo.name};
save('../brain_snRNA-seq/ConMod_SI_female_male_out.mat','filenames','Hc','modules_final','pvalues_modulePerNet','FDR2')


