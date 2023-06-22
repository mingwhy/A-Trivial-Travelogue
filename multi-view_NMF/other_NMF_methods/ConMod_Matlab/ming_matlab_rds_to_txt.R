dat=readRDS('./MF-test_co.expr/03.coexpr_net/bigscale.output/celltype_A-B-KC/net.out.rds')
dim(dat$Dp.mat)
write.table(dat$Dp.mat,file='test1.txt',quote=F,row.names = F,col.names = F)
dat=readRDS('./MF-test_co.expr/03.coexpr_net/bigscale.output/celltype_C3/net.out.rds')
dim(dat$Dp.mat)
write.table(dat$Dp.mat,file='test2.txt',quote=F,row.names = F,col.names = F)
dat=readRDS('./MF-test_co.expr/03.coexpr_net/bigscale.output/celltype_Dm8-Dm11/net.out.rds')
dim(dat$Dp.mat)
write.table(dat$Dp.mat,file='test3.txt',quote=F,row.names = F,col.names = F)
###################
dir.create('coexpr_txt')
(files=Sys.glob('../03.coexpr_net/bigscale.output/*/net.out.rds'))
for(file in files){
  cell.type=strsplit(file,'celltype_|\\/net.out.rds')[[1]][[2]]
  cat('cell.type',cell.type,'\n')
  outfile=paste0('coexpr_txt/cell.type_',cell.type,'.txt')
  dat=readRDS(file);
  write.table(dat$Dp.mat,file=outfile,quote=F,row.names = F,col.names = F)
}

###################
## in matlab
#>> A=dlmread('../test1.txt'); #remeber to add ; !!!
  
m1=dlmread('../test1.txt');
m2=dlmread('../test2.txt');
m3=dlmread('../test3.txt');

dataset = cell(3, 1);
dataset{1} = abs(m1);
dataset{2} = abs(m2);
dataset{3} = abs(m3);

% read all networks
dinfo = dir('coexpr_txt/*.txt');
dataset = cell(length(dinfo), 1);
for K = 1 : length(dinfo)
  thisfilename = dinfo(K).name;  %just the name  
  thisfilename2 = strcat('coexpr_txt/',thisfilename)
  fprintf( 'File #%d, "%s", \n', K, thisfilename2 );   %do something with the data
  dataset{K} = dlmread(thisfilename2);
end


% comment out symbol in matlab: %%
%% One-step finding conserved functional modules
tic
K = 100;
lambda = [0.01, 0.05];
xita = 2;
maxIter = 50;
num_Nodes = 2368;
modules = ConMod( dataset, num_Nodes, K, lambda, xita, maxIter );
runtime = toc;
disp(['Done.    Running time: ', num2str(runtime), ' sec.'])

# all 33, k=100, Running time: 125.0555 sec.
#Calculating the strengh matrix and the uniformity matrix...
#Obtaining candidate modules by multi-view NMF...
#K=5
#Done.    Running time: 125.0555 sec.
celldisp(modules(1,:))

%% Module validation
disp('Validation...')
permu_times = 100;
[ pvalues_modulePerNet, FDR2 ] = significantModules(modules, dataset, num_Nodes,permu_times);
disp('Done.')

out=sum(pvalues_modulePerNet,2)
size(FDR2)

filePh = fopen('output_module.txt','w');
tmp=modules{1,:};
module.id=repelem(1,size(tmp,1))
tmp(:,2)=module.id;
fprintf(filePh,'%f, %f\n',tmp'); #':t() default print by column
fclose(filePh);

filePh = fopen('output_module.txt','w');
for i = 1:size(modules,1)
  tmp=modules{i,:};
  module.id=repelem(i,size(tmp,1))
  tmp(:,2)=module.id;
  fprintf(filePh,'%f, %f\n',tmp');
end
fclose(filePh);

###################
## R code below
genes=rownames(dat$Dp.mat)
module.info=read.table('output_module.txt',sep=',')
unique(module.info[,2])

id=5
i=module.info[module.info[,2]==id,1]
genes[i]

############################################################################
## develop a R version of multiViewNMF
clear
clc
%% Load multi-network data sets
disp('Processing the networks...')

[multiNetworks, realLabels, lables_specific] = syn_dataset_overlap(0.3, false); % overlap type
num_Nodes = 500;

multiNetworks

% write matrix to files, which would be read into R
for i = 1:size(multiNetworks,1)
  A = multiNetworks{i};  
  filename =  strcat('coexpr_',string(i),'.txt');
  fid = fopen(filename,'wt');
  for ii = 1:size(A,1)
      fprintf(fid,'%g\t',A(ii,:));
      fprintf(fid,'\n');
  end
  fclose(fid)
end

[Strength, Participation] = featureNets(multiNetworks, num_Nodes);
size(Strength)
size(Participation)


% disp('Obtaining candidate modules by multi-view NMF...')
K = 5;
disp(['K=', num2str(K)])
X = {Strength, Participation};
lambda = [0.01, 0.05];
[ H, Hc, objValue ] = multiViewNMF( X, K, lambda, 50 );

X={Strength,Participation}

[multiNetworks, realLabels] = syn_dataset_common(0.5, true,'simu_dat/');

[multiNetworks, realLabels] = ming_syn_dataset_common(0.5, true,'simu_dat/');



##########################################################################
## read in Strength and Participation txt file output by R
## perform NMF inside Matlab

Strength=dlmread('../../strength.txt');
Participation=dlmread('../../participation.txt');
size(Strength)
size(Participation)

num_Nodes = size(Strength,1)

% Obtaining the candidate modules by multi-view NMF
disp('Obtaining candidate modules by multi-view NMF...')
K = 100;
disp(['K=', num2str(K)])
X = {Strength, Participation};
lambda = [0.01, 0.05];

%multiViewNMF( X, K, lambda, maxIter )
% ~1hr for 6k by 6k matrix
runtime = toc;
[ H, Hc, objValue ] = multiViewNMF( X, K, lambda, 50 ); 
disp(['Running time: ', num2str(runtime), ' sec.'])

% Selecting nodes from the consensus factors
xita = 3;
modules_final = moduleNodesSelection( Hc, xita );
runtime = toc;
disp(['Running time: ', num2str(runtime), ' sec.'])

% save modules, use R for module validation as multiNetworks is a big list of networks

%% V: Module validation
disp('Validation...') %significantModules( modules, multiNetworks, N, permu_times )
[ pvalues_modulePerNet, FDR2 ] = significantModules(modules_final, multiNetworks, num_Nodes,100);
disp('Done.')
FDR2
