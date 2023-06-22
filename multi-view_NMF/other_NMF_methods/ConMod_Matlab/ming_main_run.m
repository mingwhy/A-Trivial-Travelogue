
clear
clc

% I: Automatically generated synthetic datasets
% create two folder first 'simu_dat_common' and 'simu_dat_overlap'
[multiNetworks, realLabels] = ming_syn_dataset_common(0.5, true,'simu_dat_common/'); %50nets with 500 nodes each
[multiNetworks, realLabels, lables_specific] = ming_syn_dataset_overlap(0.3, true,'simu_dat_overlap/'); % 15 nets with 500 noedes each, overlap type


%% II: Load multi-network data sets
disp('Processing the networks...')

% Load data from files
nets = importdata('./simu_dat_overlap/networklist.txt');
T = length(nets);
multiNetworks = cell(T, 1);
for i = 1:T
    disp(['network: ', num2str(i)])
    multiNetworks{i} = load(['./simu_dat_overlap/', nets{i}]);
end
xrealLabels = importdata('./simu_dat_overlap/labels.txt');
for(i =1:size(realLabels,1))
    realLabels{i}=xrealLabels(i,:);
end
size(multiNetworks{1})


%% III: One-step finding conserved functional modules
tic
K = 5;
lambda = [0.01, 0.05];
xita = 2;
maxIter = 50;
num_Nodes = 500;
modules = ConMod( multiNetworks, num_Nodes, K, lambda, xita, maxIter );
runtime = toc;
disp(['Done.    Running time: ', num2str(runtime), ' sec.'])

%% Clustering performance
[ TPR, FPR, Accuracy, MCC] = evaluation(modules, realLabels, num_Nodes);
TPR
FPR 
Accuracy
MCC


%% Alternatively, Step-by-step finding conserved functional modules
% Calculting the feature networks
tic
disp('Calculating the strengh matrix and the uniformity matrix...')
num_Nodes = 500;
[Strength, Participation] = featureNets(multiNetworks, num_Nodes);
 
% Obtaining the candidate modules by multi-view NMF
disp('Obtaining candidate modules by multi-view NMF...')
K = 5;
disp(['K=', num2str(K)])
X = {Strength, Participation};
lambda = [0.01, 0.05];
[ H, Hc, objValue ] = multiViewNMF( X, K, lambda, 50 ); %multiViewNMF( X, K, lambda, maxIter )
 
% Selecting nodes from the consensus factors
xita = 1.5;
modules_final = moduleNodesSelection( Hc, xita );
runtime = toc;
disp(['Running time: ', num2str(runtime), ' sec.'])


%% IV: Clustering performance
[ TPR, FPR, Accuracy, MCC] = evaluation(modules_final, realLabels, num_Nodes);



%% V: Module validation
disp('Validation...') %significantModules( modules, multiNetworks, N, permu_times )
[ pvalues_modulePerNet, FDR2 ] = significantModules(modules_final, multiNetworks, num_Nodes,100);
disp('Done.')
FDR2

% for 'common' type simulated data, number of 0 in each row
% C1 exist in 5 out of 30 simulated matrix
% C2 exist in 25
% C3 exist in 20
% C4 exist in 10
% C5 exist in 15
