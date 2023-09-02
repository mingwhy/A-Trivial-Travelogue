% flux sampling: https://github.com/HAHerrmann/FluxSamplingComparison
%%%%%%%%
%load model
load('/Users/mingyang/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/Fruitfly-GEM.mat');
model=fruitflyGEM;

load('body_lb_ub.mat') 
size(imat_ub) 

%% Sampling CHRR
for i=[1,2]
   
    model1 = model;
    model1.ub = imat_ub(:,1);
    model1.lb =  imat_lb(:,1);
    
    %https://github.com/ankueken/Chlamy_model/blob/master/sampling_sst_f.m
    %options.nFiles = 5; %ACHR only.
    options.nPointsPerFile = 100;
    options.nPointsReturned = 100;
    
    %set constraint bounds
    model1.b = zeros(size(model1.mets));
   
   %https://opencobra.github.io/cobratoolbox/latest/modules/analysis/sampling/index.html
   %[modelSampling, samples] = sampleCbModel(model, sampleFile, samplerName, options, modelSampling)
   
   tic
   [Model_reduced, samples] =  sampleCbModel(model1, [],'CHRR', options); %outputs model (polytope) and samples 
   toc
   
   filenameA = insertBefore('ResultsCHRR__A_T_2.mat',17,num2str(sDens));
   filename = insertBefore(filenameA,13,modelname);
   filenameB = insertBefore('ResultsCHRR__A_T_12.csv',17,num2str(sDens));
   filename1 = insertBefore(filenameB,13,modelname);
   %save(filename, 'samples')
   t_samples = transpose(samples);
   csvwrite(filename1,t_samples)
   disp('Finished model')
