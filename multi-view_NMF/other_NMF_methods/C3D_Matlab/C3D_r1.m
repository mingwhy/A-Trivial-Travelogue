function C3D_r1(file_path,qval,vec_n,vad_tag,varargin)
% C3D: main function
%
% Input (mandatory):
% The () after the variable names indicate how the variable is
% called in the instructions (see details in "Instructions for C3D.pdf") 
% file_path (directory): specifies the location of the mutiple input data sets 
% qval (MER): threshold to be used for cluster nodes selection (a float number in [0,1])
% vec_n (n.vectors): the number of candidate vectors
% vad_tag (sim.indicator): a switch for cluster validation
%
% Input (optional):
% varargin
%
% Written by Xiaolin Xiao, 2012


pwd; 
sourcefolder=pwd; 


sf=dir('vx.txt');
if size(sf,1)==1
    delete('vx.txt');
end
sf=dir('vx_outqval_all.txt');
if size(sf,1)==1
    delete('vx_outqval_all.txt');
end

c=clock; 
dfile=['run_c3d_main_at_y' int2str(c(1)) '_m' int2str(c(2)) '_d' int2str(c(3)) '_h' int2str(c(4)) '_min' int2str(c(5)) '_sec' int2str(c(6)) '.txt'];
cd(file_path);
mkdir('results');
cd('results');
diary(dfile);
diary on;
cd(sourcefolder);
fprintf('start running: year%d month%d day%d hour%d min%d sec%d\n',c(1),c(2),c(3),c(4),c(5),c(6));

% input parameters check 
if nargin<1
    error(message('TooFewInputs'));
end
if (nargout > 8)
    error(message('MATLAB:clustervalidation:TooManyOutputs'));
end

if ~isnumeric(qval) || ~isreal(qval)
    error('the threshold for MER must be a real number');
end
if qval>1 || qval<0 
    error('the threshold for MER must be a non-negative real number lower or equal to 1');
end

if ~isnumeric(vec_n) || ~isreal(vec_n)
    error('the number of candidate vectors MUST be a real number');
end
if vec_n<0 || mod(vec_n,1)~=0
    error('the number of candidate vectors MUST be a positive interger');
end

% load data sets
tic;
cd(file_path); 
file='datasets_list.txt'; 
fid=fopen(file,'r');
cond=0;
while ~feof(fid)
    fscanf(fid,'%s\n',1);
    cond=cond+1;
end
fclose(fid);
if cond<2
    error('number of conditions is less than 2 or mistake in file datasets_list.txt');
end
fprintf('Step 1: loading %d data sets......\n',cond);
D=cell(cond,1);
fid=fopen(file,'r');
vers=version('-release');
vers=vers(1:4);
vers=str2double(vers);
for i=1:cond
    a=fscanf(fid,'%s\n',1);
    filename=a; 
    disp(filename);
    if vers>2012
        D{i}=load(filename);
    else
        D{i}=importdata(filename);
    end
end
fclose(fid);

file='geneid.txt';
fid=fopen(file,'r');
N=0;
while ~feof(fid)
    fscanf(fid,'%s\n',1);
    N=N+1;
end
fclose(fid);
if N<=0
    error('mistake: node names in file geneid.txt');
end
genes=cell(N,1);
fid=fopen(file,'r');
for i=1:N
    genes{i}=fscanf(fid,'%s\n',1);
end
fclose(fid);
time_int=toc/60;
fprintf('Step 1 Done: %f minutes\n',time_int);


cd('results');
cdqval=['fdr' num2str(qval) 'cut'];
mkdir(cdqval);
sp_cdqval=[cdqval '/clusters'];
mkdir(sp_cdqval);
cd(cdqval);
mkdir('edgelist');    

cd(sourcefolder); 


if cond~=length(D)
    error('Data size mismatch');
end
drsize=nan(1,cond);
real_tag=nan(1,cond);
ep_tag=nan(1,cond);
for i=1:cond
    drsize(i)=size(D{i},1);
    real_tag(i)=isreal(D{i});
    ep_tag(i)=isempty(D{i});
    if sum(sum(isnan(D{i})))>0
       error('nan in data');
    end
    if sum(sum(isinf(D{i})))>0
        error('Inf in data');
    end    
end
if ~isequal(drsize, drsize(1).*ones(1,cond))
    error(message('MATLAB:hogsvd:MatrixRowMismatch')) 
end
if ~isequal(drsize, N.*ones(1,cond))
    error('The %d nodes your provided are NOT in ROWS, check your data (tranpose the data if necessary) and re-try\n',N);
end
if ~isequal(real_tag,ones(1,cond))
    error('Complex number exist');
end
if ~isequal(ep_tag,zeros(1,cond))
    error('Empty data matrix exist');
end

issym = @(x) isequal(x,x.'); 
tf=nan(cond,1);
for i=1:cond
    tf(i)=issym(D{i});
end
if isequal(tf,ones(cond,1)) 
    data_type=1; 
else
    data_type=0;
end

% display optional varagrin
[std,crtag,rnd,rnde,thrvp,ctrig,rab2,abst]=checkInputs(varargin{:});
disp('Here are the optional input arguments for preprocessing the data:');
fprintf('standardization method=%d,correlation tag=%d\n',std,crtag);
% data preprocessing
if std>0
    disp('Data preprocessing (Standardization using zscore)......');
    if std==1 
        for i=1:cond
            D{i}=zscore(D{i}');
            D{i}=D{i}';
            if sum(sum(isnan(D{i})))>0
               error('nan in data (zscored), adjust your data or re-try with different arguments');
            end
            if sum(sum(isinf(D{i})))>0
                error('Inf in data (zscored), adjust your data or re-try with different arguments');
            end
            if ~isreal(D{i})
                error('Complex number in data (zscored), adjust your data or re-try with different arguments');
            end
        end       
    end
    if std==2 
        for i=1:cond
            D{i}(D{i}==0)=1; 
            D{i}=log2(D{i}); 
            if sum(sum(isnan(D{i})))>0
               error('nan in data (log scaled), adjust your data or re-try with different arguments');
            end
            if sum(sum(isinf(D{i})))>0
                error('Inf in data (log scaled), adjust your data or re-try with different arguments');
            end 
            if ~isreal(D{i})
                error('Complex number in data (log scaled), adjust your data or re-try with different arguments');
            end
        end
    end        
end
% correlation 
if crtag>0
    disp('Data preprocessing (Computing correlation matrices)......');
    data_type=1; 
    if crtag==1 
        disp('Generating Pearson Correlation for the larg size data may take a lot of time...');
        for i=1:cond
            D{i}=corrcoef(D{i}');
            if sum(sum(isnan(D{i})))>0
               error('nan in data (Pearson Correlation), adjust your data or re-try with different arguments');
            end
            if sum(sum(isinf(D{i})))>0
                error('Inf in data (Pearson Correlation), adjust your data or re-try with different arguments');
            end 
            if ~isreal(D{i})
                error('Complex number in data (Pearson Correlation), adjust your data or re-try with different arguments');
            end
        end   
    end
    if crtag==2 
        disp('Generating Kendall Correlation for your data may take a lot of time...'); 
        tic;
        for i=1:cond
            D{i}=corr(D{i}','type','Kendall'); 
            if sum(sum(isnan(D{i})))>0
               error('nan in data (Kendall Correlation), adjust your data or re-try with different arguments');
            end
            if sum(sum(isinf(D{i})))>0
                error('Inf in data (Kendall Correlation), adjust your data or re-try with different arguments');
            end 
            if ~isreal(D{i})
                error('Complex number in data (Kendall Correlation), adjust your data or re-try with different arguments');
            end
        end
        time_int=toc/60;
        fprintf('Computing Kendall Correlation Done: %f minutes\n',time_int);
    end
    if crtag==3 
        disp('Generating Spearman Correlation for your data may take a lot of time...');
        tic;
        for i=1:cond
            D{i}=corr(D{i}','type','Spearman');
            if sum(sum(isnan(D{i})))>0
               error('nan in data (Spearman Correlation), adjust your data or re-try with different arguments');
            end
            if sum(sum(isinf(D{i})))>0
                error('Inf in data (Spearman Correlation), adjust your data or re-try with different arguments');
            end 
            if ~isreal(D{i})
                error('Complex number in data (Spearman Correlation), adjust your data or re-try with different arguments');
            end
        end 
        time_int=toc/60;
        fprintf('Computing Spearman Correlation Done: %f minutes\n',time_int);
    end    
end

tic;
disp('Step 2: Compute the HO GSVD decomposition......');
if cond==2
    Sigma=cell(2,1);
    [~,~,V,Sigma{1},Sigma{2}]=gsvd(D{1}',D{2}',0); 

    if isequal(size(Sigma{1}),size(Sigma{2})) && isequal(size(Sigma{1},1),size(Sigma{1},2)) 
        Lambda=diag(diag(Sigma{1})./diag(Sigma{2})); 
    else
        lam1=zeros(size(V,2),1); 
        lam2=zeros(size(V,2),1);
        mins1=min(size(Sigma{1},1),size(Sigma{1},2));
        lam1(end-mins1+1:end)=diag(Sigma{1}(end-mins1+1:end-mins1+1,end-mins1+1:end-mins1+1));
        mins2=min(size(Sigma{2},1),size(Sigma{2},2));
        lam2(1:mins2)=diag(Sigma{2}(1:mins2,1:mins2));
        lam2(mins2+1:end)=lam2(mins2)./(1e16); 
        Lambda=diag(lam1./lam2);     
    end
    
end

if cond>2
    if data_type==1 % for co-expression data
        [V,Lambda]=e_trans(D); 
    else
        [V,Lambda]=hogsvd_trans(D); 
    end
end
time_int=toc/60;
fprintf('Step 2 Done: %f minutes\n',time_int);

if vec_n>size(V,2)
    fprintf('The given number of candidate vectors (%d) exceeds column size of V, here change to use all %d vectors in V!\n',vec_n,size(V,2));
    vec_n=size(V,2);
end

lam=diag(Lambda);
if sum(sum(isnan(V)))>0 || sum(sum(isinf(V)))>0
    error('Nan or Inf in V');
end
if sum(isnan(lam))>0 || sum(isinf(lam))>0    
    error('Nan or Inf in the diagonal entries');
end


dist=lam-1;
dist=abs(dist);
[~,key1]=sort(dist); 
[~,eigsdx]=sort(abs(lam)); 

if cond>2 
    if size(V,2)>=vec_n
        key=key1(1:vec_n);
    else 
        key=key1;
    end
else
    if cond==2 
        if size(V,2)>=vec_n 
            if round(vec_n/3)>5
                key=union(eigsdx(1:5),eigsdx(end-4:end)); 
                if length(key)>10 || length(key)>=vec_n
                    error('Mistake: length of candidate vectors');
                else
                    key=union(key1(1:vec_n-length(key)),key); 
                end               
            else
                key=union(eigsdx(1:round(vec_n/3)), eigsdx(end-round(vec_n/3)+1:end));
                if length(key)>2*round(vec_n/3) || length(key)>=vec_n
                    error('Mistake: length of candidate vctors');
                else 
                    key=union(key1(1:vec_n-length(key)),key); 
                end               
            end           
        else
            key=key1;
        end
    else
        error('Mistake: size of V');
    end
end


if length(key)<vec_n 
    [~,dif_pos]=setdiff(key1,key); 
    for i=1:vec_n-length(key)
        key=union(key,key1(dif_pos(i))); 
    end   
end


if length(key)~=vec_n
    error('Mismatch: size of selected vectors');
end

fprintf('The following %d candidate vectors are selected:\n',length(key));
disp(key(1:end)); 


if ~isempty(key) 
    tic;
    disp('Step 3: Calling fdrtool to select the cluster nodes......');

    cd(sourcefolder);
    
    filen=0;
    c_name=cell(2*length(key),1); 
    v_num=nan(2*length(key),1); 
    overlap_m=nan(length(c_name));
    overlap_node=cell(length(c_name));
    
    cd(file_path);
    cd('results'); 
    cd(cdqval);
    file2='Number_all_clusters_overlap.txt'; 
    fid2=fopen(file2,'wt');
    file3='Nodes_all_clusters_overlap.txt'; 
    fid3=fopen(file3,'wt');
    
    
    for k=1:vec_n
        cd(sourcefolder);
        vx=V(:,key(k));
        save('vx.txt','vx','-ascii','-double','-tabs');
        disp('**************************************************************');
        fprintf('Selecting cluster nodes by fdrtool for v%d\n',key(k));
        !unset DYLD_LIBRARY_PATH; /usr/bin/Rscript fdrqval.r
      
        
        sf=dir('vx_outqval_all.txt');
        if size(sf,1)==1 
            if sf.bytes~=0
                cqval=importdata('vx_outqval_all.txt');       
                if isfield(cqval,'data')
                    cqval=cqval.data;
                    qval_all=cqval; 
               
                    cluster_pos=1:1:N; 
                    cluster=genes(cluster_pos); 
       
                    cd(file_path); 
                    cd('results'); 
                    cd(cdqval);
                    cluster_pos=cluster_pos(cqval<=qval); 
                    cluster=cluster(cqval<=qval);
                    cqval=cqval(cqval<=qval);
                    fprintf('cut %d nodes for this cluster (MER=%f)\n',length(cluster),qval);
                    if ~isequal(length(cluster),length(cqval))
                        error('Size mismatch: the cluster selected by MER=%f\n',qval);
                    end
                    fout=['v' int2str(key(k)) '_c' int2str(length(cluster)) 'nodes_q' num2str(qval) '.txt'];
                    fid=fopen(fout,'wt');
                    if isequal(length(cluster),length(cqval))
                        for m=1: length(cluster)
                            fprintf(fid,'%s\t%f\t%d\n',cluster{m},cqval(m),cluster_pos(m)); 
                        end
                        fclose(fid);
                    else
                        error('Mismatch: length of current output cluster and the corresponding MER-value');
                    end

                    if length(cluster)<N && ~isempty(cluster)
                        [~,b]=sort(V(:,key(k)));
                        [~,~,ib]=intersect(cluster_pos,b);
                        if length(ib)==N 
                            fprintf('Warning: using v%d,the cluster selected covers all nodes ... re-try with a different MER or adjust your data\n',key(k));
                        end        

                        cd('clusters'); 
                        ful=1:1:N;
                        ac=setdiff(ful,sort(ib)); 
                        if ac(1)==1
                            gene_pos_b=b(ac(end)+1:end); 
                            clusterbot=genes(gene_pos_b);
                            clusterbot_q=qval_all(b(ac(end)+1:end));
                            filenameb=['v' int2str(key(k)) '_bot' int2str(length(clusterbot)) '_c' int2str(length(cluster)) 'nodes_q' num2str(qval) '.txt']; 
                            fid=fopen(filenameb,'wt');
                            for m=1:length(clusterbot)
                                fprintf(fid,'%s\t%f\t%d\n',clusterbot{m},clusterbot_q(m),gene_pos_b(m));
                            end 
                            fclose(fid);

                           
                            filen=filen+1;
                            c_name{filen}=filenameb;
                            v_num(filen)=key(k);
                            cd(file_path);
                            cd('results');
                            cd(cdqval); 
                            fprintf(fid2,'%s\t',filenameb);
                            fprintf(fid3,'%s\t',filenameb);

                            
                            cd(file_path);
                            cd('results');
                            cd(cdqval);
                            cd('edgelist');
                            corv=nan(cond,1);
                            file_edg=['edge_v' int2str(key(k)) '_bot' int2str(length(clusterbot)) '_c' int2str(length(cluster)) 'nodes_q' num2str(qval) '.txt'];
                            fid=fopen(file_edg,'wt');
                            for i=1:length(clusterbot)
                                for j=i+1:length(clusterbot)
                                    fprintf(fid,'%s\t%s\t',clusterbot{i},clusterbot{j});
                                    if crtag==0 
                                        for n=1:cond
                                            corv(n)=sum(D{n}(gene_pos_b(i),:).*D{n}(gene_pos_b(j),:)); 
                                            fprintf(fid,'%f\t',corv(n));
                                        end
                                        fprintf(fid,'\n');
                                    else
                                        for n=1:cond
                                            corv(n)=D{n}(i,j);
                                            fprintf(fid,'%f\t',corv(n));
                                        end
                                        fprintf(fid,'\n');
                                    end
                                end
                            end
                            fclose(fid);

                        else
                            if ac(end)==N
                                gene_pos_t=b(1:ac(1)-1);
                                clustertop=genes(gene_pos_t);
                                clustertop_q=qval_all(b(1:ac(1)-1));
                                filenamet=['v' int2str(key(k)) '_top' int2str(length(clustertop)) '_c' int2str(length(cluster)) 'nodes_q' num2str(qval) '.txt'];
                                fid=fopen(filenamet,'wt');
                                for m=1:length(clustertop)
                                    fprintf(fid,'%s\t%f\t%d\n',clustertop{m},clustertop_q(m),gene_pos_t(m));
                                end
                                fclose(fid);

                                
                                filen=filen+1;
                                c_name{filen}=filenamet;
                                v_num(filen)=key(k);
                                cd(file_path);
                                cd('results');
                                cd(cdqval);
                                fprintf(fid2,'%s\t',filenamet);
                                fprintf(fid3,'%s\t',filenamet);

                                cd(file_path);
                                cd('results');
                                cd(cdqval);
                                cd('edgelist');
                                corv=nan(cond,1);
                                file_edg=['edge_v' int2str(key(k)) '_top' int2str(length(clustertop)) '_c' int2str(length(cluster)) 'nodes_q' num2str(qval) '.txt'];
                                fid=fopen(file_edg,'wt');
                                for i=1:length(clustertop)
                                    for j=i+1:length(clustertop) 
                                        fprintf(fid,'%s\t%s\t',clustertop{i},clustertop{j});
                                        if crtag==0
                                            for n=1:cond
                                                corv(n)=sum(D{n}(gene_pos_t(i),:).*D{n}(gene_pos_t(j),:)); 
                                                fprintf(fid,'%f\t',corv(n));
                                            end
                                            fprintf(fid,'\n');
                                        else 
                                            for n=1:cond
                                                corv(n)=D{n}(i,j); 
                                                fprintf(fid,'%f\t',corv(n));
                                            end
                                            fprintf(fid,'\n');
                                        end
                                    end
                                end
                                fclose(fid);

                            else
                                gene_pos_t=b(1:ac(1)-1);
                                clustertop=genes(gene_pos_t);
                                clustertop_q=qval_all(b(1:ac(1)-1));
                                filenamet=['v' int2str(key(k)) '_top' int2str(length(clustertop)) '_c' int2str(length(cluster)) 'nodes_q' num2str(qval) '.txt'];
                                fid=fopen(filenamet,'wt');
                                for m=1:length(clustertop)
                                    fprintf(fid,'%s\t%f\t%d\n',clustertop{m},clustertop_q(m),gene_pos_t(m));
                                end
                                fclose(fid);

                               
                                filen=filen+1;
                                c_name{filen}=filenamet;
                                v_num(filen)=key(k);
                                rc=filen;
                                cd(file_path);
                                cd('results');
                                cd(cdqval);
                                fprintf(fid2,'%s\t',filenamet);
                                fprintf(fid3,'%s\t',filenamet);

                                
                                cd('clusters');
                                gene_pos_b=b(ac(end)+1:end);
                                clusterbot=genes(gene_pos_b);
                                clusterbot_q=qval_all(b(ac(end)+1:end));
                                filenameb=['v' int2str(key(k)) '_bot' int2str(length(clusterbot)) '_c' int2str(length(cluster)) 'nodes_q' num2str(qval) '.txt'];
                                fid=fopen(filenameb,'wt');
                                for m=1:length(clusterbot)
                                    fprintf(fid,'%s\t%f\t%d\n',clusterbot{m},clusterbot_q(m),gene_pos_b(m)); 
                                end
                                fclose(fid);

                                
                                filen=filen+1;
                                c_name{filen}=filenameb;
                                v_num(filen)=key(k);
                                sc=filen;
                                cd(file_path);
                                cd('results');
                                cd(cdqval);
                                fprintf(fid2,'%s\t',filenameb);
                                fprintf(fid3,'%s\t',filenameb);

                                overlap_m(rc,sc)=0; 
                                overlap_m(sc,rc)=0; 

                                cd(file_path);
                                cd('results');
                                cd(cdqval);
                                cd('edgelist');
                                file_edg=['edge_v' int2str(key(k)) '_top' int2str(length(clustertop)) '_c' int2str(length(cluster)) 'nodes_q' num2str(qval) '.txt'];
                                fid=fopen(file_edg,'wt');
                                for i=1:length(clustertop)
                                    for j=i+1:length(clustertop)
                                        fprintf(fid,'%s\t%s\t',clustertop{i},clustertop{j});
                                        if crtag==0
                                            for n=1:cond
                                                corv(n)=sum(D{n}(gene_pos_t(i),:).*D{n}(gene_pos_t(j),:));
                                                fprintf(fid,'%f\t',corv(n));
                                            end
                                            fprintf(fid,'\n');
                                        else 
                                            for n=1:cond
                                                corv(n)=D{n}(i,j);
                                                fprintf(fid,'%f\t',corv(n));
                                            end
                                            fprintf(fid,'\n');
                                        end
                                    end
                                end
                                fclose(fid);

                                corv=nan(cond,1);
                                file_edg=['edge_v' int2str(key(k)) '_bot' int2str(length(clusterbot)) '_c' int2str(length(cluster)) 'nodes_q' num2str(qval) '.txt'];
                                fid=fopen(file_edg,'wt');
                                for i=1:length(clusterbot)
                                    for j=i+1:length(clusterbot)
                                        fprintf(fid,'%s\t%s\t',clusterbot{i},clusterbot{j});
                                        if crtag==0
                                            for n=1:cond
                                                corv(n)=sum(D{n}(gene_pos_b(i),:).*D{n}(gene_pos_b(j),:));
                                                fprintf(fid,'%f\t',corv(n));
                                            end
                                            fprintf(fid,'\n');
                                        else 
                                            for n=1:cond
                                                corv(n)=D{n}(i,j); 
                                                fprintf(fid,'%f\t',corv(n));
                                            end
                                            fprintf(fid,'\n');
                                        end
                                    end
                                end
                                fclose(fid);

                            end
                        end
                    else
                        disp('An empty cluster or the cluster covers all nodes: no cluster/edge list generated, no validation');  
                    end 
                        clear cluster_pos;       
                        clear cluster;
                                    
                else
                    fprintf('v%d: no real data (empty output) in output files generated from fdrtool\n',key(k));
                end 
                
                clear cqval;
                
            else 
                fprintf('v%d: empty output files generated from fdrtool\n',key(k)); 
            end
            
            cd(sourcefolder);
            delete('vx.txt'); 
            delete('vx_outqval_all.txt'); 
            
        else
            delete('vx.txt'); 
            fprintf('v%d: no output files generated from fdrtool\n',key(k));
        end
    end 
    
    
    cd(file_path);
    cd('results');
    cd(cdqval);
    fprintf(fid2,'\n');
    fprintf(fid3,'\n');
    
    if filen<=length(overlap_m)
        overlap_m=overlap_m(1:filen,1:filen);
        overlap_node=overlap_node(1:filen,1:filen);
    else
        error('overlap matrix size mismatch');
    end
    
    
    v_num(filen+1:end)=[];
    if length(v_num)~=filen
        error('Mistake: total number of cluster files checked');
    end
    
    
    for r=1:filen
        for s=r+1:filen
            % file r
            filer=c_name{r};
            headr_len=1+length(int2str(v_num(r)))+4;  
            % file s
            files=c_name{s};
            heads_len=1+length(int2str(v_num(s)))+4; 

            if isnan(overlap_m(r,s))
                
                tailr=headr_len+1;
                while ~strcmp(filer(tailr+1),'_') 
                    tailr=tailr+1;
                end
                clen=[filer(headr_len+1)];
                if tailr>headr_len+1 
                    for i=headr_len+2:tailr
                        clen=[clen filer(i)]; 
                    end
                end
                clen=str2double(clen);
                cd(file_path);
                cd('results');
                cd(cdqval);
                cd('clusters');
                fid=fopen(filer,'r');
                rc=cell(clen,1);
                for i=1:clen
                    rc{i}=fscanf(fid,'%s\t',1); 
                    fscanf(fid,'%s\n',2); 
                end
                fclose(fid);

                % file s          
                tails=heads_len+1;
                while ~strcmp(files(tails+1),'_') 
                    tails=tails+1;
                end
                clen=[files(heads_len+1)];
                if tails>heads_len+1 
                    for i=heads_len+2:tails
                        clen=[clen files(i)];
                    end
                end
                clen=str2double(clen);
                cd(file_path);
                cd('results');
                cd(cdqval);
                cd('clusters');
                fid=fopen(files,'r');
                sc=cell(clen,1);
                for i=1:clen
                    sc{i}=fscanf(fid,'%s\t',1); 
                    fscanf(fid,'%s\n',2);
                end
                fclose(fid);

                % compare file r and file s (2 clusters)
                overlap_node{r,s}=intersect(rc,sc);  
                overlap_m(r,s)=length(overlap_node{r,s}); 
            end
        end
    end
    
    overlap_m=triu(overlap_m,1);
    overlap_m=overlap_m+overlap_m';     
        
    
    for j=1:filen
        filename=c_name{j};
        head_len=1+length(int2str(v_num(j)))+4; 
        tail=head_len+1;
        while ~strcmp(filename(tail+1),'_') 
            tail=tail+1;
        end
        clen=[filename(head_len+1)];
        if tail>head_len+1 
            for i=head_len+2:tail
                clen=[clen filename(i)]; 
            end
        end
        overlap_m(j,j)=str2double(clen);
        overlap_node{j,j}='all';
    end 
    
    cd(file_path);
    cd('results');
    cd(cdqval);
    % write overlap
    for r=1:filen
        for s=1:filen
            fprintf(fid2,'%d\t',overlap_m(r,s));
            if r~=s
                tt=overlap_node{r,s};
            else
                tt=overlap_node(r,s); 
            end
            for j=1:length(tt)
                fprintf(fid3,'%s ',tt{j}); 
            end
            fprintf(fid3,'\t'); 
        end
        fprintf(fid2,'\n');
        fprintf(fid3,'\n');
    end
    
    fclose(fid2);
    fclose(fid3);
    
    time_int=toc/60;
    fprintf('Step 3 Done: %f minutes\n',time_int);
    
else
    error('Mistake: empty candidate vectors');
end


cd(file_path);
cd('results');
cd(cdqval);
delete('v*');

if vad_tag==1 

    cd(sourcefolder);
    tic;
    disp('Step 4: Validating all the clusters detected above......');
    
    disp('Here are the arguments used for cluster validation:');
    fprintf('randomization times=%d to %d, threshold=%f, individual cluster quality measure=%d, overall cluster quality measure=%d\n',rnde,rnd,thrvp,ctrig,rab2);
    if abst==1;
        disp('And use the absolute values from the data in cluster validation');
    end

    cd(file_path);
    cd('results'); 
    if abst==0
        file2=['clusters_summary_q' num2str(qval) '_swap' int2str(rnde) 'to' int2str(rnd) '_thrv' int2str(thrvp*100) '_percq' int2str(ctrig) '_ocq' int2str(rab2) '.txt']; 
    else 
        file2=['clusters_summary_q' num2str(qval) 'ABS_swap' int2str(rnde) 'to' int2str(rnd) '_thrv' int2str(thrvp*100) '_percq' int2str(ctrig) '_ocq' int2str(rab2) '.txt']; 
    end   
    
    fid2=fopen(file2,'wt');
    fprintf(fid2,'cluster\toverall p\tconditions\toverall randomization\toverall cluster quality\t');
    for i=1:cond
        fprintf(fid2,'individual p in condition%d\t',i);
        fprintf(fid2,'individual randomization in condition%d\t',i);
        fprintf(fid2,'individual cluster quality in condition%d\t',i);
    end
    fprintf(fid2,'cluster size\n');
    
    for k=1:vec_n 
        fprintf('validating clusters selected from v%d\n',key(k));
        cd(file_path); 
        cd('results');
        cd(cdqval);
        cd('clusters');
        
        
        filenamet=['v' int2str(key(k)) '_top*'];
        sf=dir(filenamet);
        if size(sf,1)==1 
            if sf.bytes~=0 
                filenamet=sf.name;
                disp('**************************************************************');
                fprintf('now validating the cluster (file) %s\n',filenamet);
                disp('******************************************');
                fid=fopen(filenamet,'r');
                clen=0;
                while ~feof(fid)
                    fscanf(fid,'%s\n',3); 
                    clen=clen+1;
                end
                fclose(fid);
                td1=1;
                td2=clen;   
                if clen>=3 && clen<N 
                    cd(sourcefolder); 
                    [ds,pv,itr_flag,rab,pvrab,itr_flagrab,consout]=validate_trans(data_type,D,V,key(k),td1,td2,rnd,rnde,thrvp,ctrig,rab2,abst);
                    disp('******************************************');
                    time_int=toc/60;
                    fprintf('Validation for cluster (%s) is done: %f minutes\n',sf.name,time_int);

                    cd(file_path);
                    cd('results');
                    
                    fprintf(fid2,'%s\t%f\t%s\t%d\t%f\t',filenamet,pvrab,int2str(consout),itr_flagrab,rab);
                    for i=1:cond
                        fprintf(fid2,'%f\t',pv(i));
                        fprintf(fid2,'%d\t',itr_flag(i));
                        fprintf(fid2,'%f\t',ds(i));
                    end
                    fprintf(fid2,'%d\n',td2-td1+1);
                else
                    delete(filenamet);
                    disp('This cluster is too small(<3 nodes)--delete');
                end
   
            end
        end
        
        cd(file_path); 
        cd('results');
        cd(cdqval);
        cd('clusters');
        filenameb=['v' int2str(key(k)) '_bot*'];
        sf=dir(filenameb);
        if size(sf,1)==1 
            if sf.bytes~=0 
                
                filenameb=sf.name;
                disp('**************************************************************');
                fprintf('now validating the cluster (file) %s\n',filenameb);
                disp('******************************************');
                
                fid=fopen(filenameb,'r');
                clen=0;
                while ~feof(fid)
                    fscanf(fid,'%s\n',3);
                    clen=clen+1;
                end
                fclose(fid);
                td1=N-clen+1;
                td2=N; 
                if clen>=3 && clen<N
                    cd(sourcefolder);                
                    [ds,pv,itr_flag,rab,pvrab,itr_flagrab,consout]=validate_trans(data_type,D,V,key(k),td1,td2,rnd,rnde,thrvp,ctrig,rab2,abst);              
                    disp('******************************************');
                    time_int=toc/60;
                    fprintf('Validation for cluster (%s) is done: %f minutes\n',sf.name,time_int);

                    cd(file_path);
                    cd('results');
                    
                    fprintf(fid2,'%s\t%f\t%s\t%d\t%f\t',filenameb,pvrab,int2str(consout),itr_flagrab,rab);
                    for i=1:cond
                        fprintf(fid2,'%f\t',pv(i));
                        fprintf(fid2,'%d\t',itr_flag(i));
                        fprintf(fid2,'%f\t',ds(i));
                    end
                    fprintf(fid2,'%d\n',td2-td1+1);
                else
                    delete(filenameb);
                    disp('This cluster is too small(<3 nodes)--delete');
                    disp('******************************************');
                end
           
            end
        end
        
    end     
    
    fclose(fid2);

    time_int=toc/60;
    fprintf('Step 4 Done: %f minutes\n',time_int);
end

diary off;



%-------------------------------------------------------------------------%
% Nested functions
%-------------------------------------------------------------------------%
function [std,crtag,rnd,rnde,thrvp,ctrig,rab2,abst] = checkInputs(varargin) 
% Nested funtion for optional input varaibles
% Instructions:
% std      = normalize
% crtag    = correlation
% max.perm = rnd
% min.perm = rnde
% thrvp    = p.threshold
% ctrig    = individual.cluster.quality
% rab2     = overall.cluster.quality
% abst     = norm
%
% Written by Xiaolin Xiao, 2012


 if nargin>8
    error('Too many input arguments');
 end
 
 std=0;
 crtag=0;
 rnd=1000;
 rnde=100;
 thrvp=0.05;
 ctrig=2;
 rab2=5;
 abst=0;
 
 if ~isempty(varargin) 
    switch varargin{1}
        case 'zscore'
            std=1;
        case 'log'
            std=2;
        otherwise
            switch varargin{1}
                case 'Pearson'
                   crtag=1;
                case 'Kendall'
                   crtag=2; 
                case 'Spearman'
                   crtag=3;
                otherwise
                   if isfloat(varargin{1}) && mod(varargin{1},1)==0 && varargin{1}>=100 
                       if length(varargin)>1
                           if isfloat(varargin{2}) && mod(varargin{2},1)==0 && varargin{2}>=100
                               rnd=max(varargin{1},varargin{2}); 
                               rnde=min(varargin{1},varargin{2}); 
                               if length(varargin)>2
                                   if isfloat(varargin{3})
                                       if varargin{3}>0 && varargin{3}<1 
                                           thrvp=varargin{3};
                                           if length(varargin)>3    
                                               if ischar(varargin{4})
                                                   if varargin{4}(1)=='c'
                                                       if strcmp(varargin{4},'c1') || strcmp(varargin{4},'c2')
                                                           ctrig=str2double(varargin{4}(2));
                                                           if length(varargin)>4
                                                               if varargin{5}(1)=='q' && mod(str2double(varargin{5}(2)),1)==0
                                                                   if str2double(varargin{5}(2))<=7 && str2double(varargin{5}(2))>=1 && length(varargin{5})<3 
                                                                       rab2=str2double(varargin{5}(2)); 
                                                                       if length(varargin)>5
                                                                           if strcmp(varargin{6},'abs') 
                                                                               abst=1;
                                                                               if length(varargin)>6
                                                                                   disp('input after the end string (abs) will be ignored');
                                                                               end
                                                                           else
                                                                               error('unidentified input characters or improper order for the optional input arguments');
                                                                           end  
                                                                       end
                                                                   else
                                                                       error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                                   end
                                                               else
                                                                   if length(varargin)>4
                                                                       if strcmp(varargin{5},'abs')
                                                                           abst=1;
                                                                           if length(varargin)>5
                                                                               disp('input after the end string (abs) will be ignored');
                                                                           end
                                                                       else
                                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                                       end
                                                                   end
                                                               end
                                                           end        
                                                       else
                                                           error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                                       end
                                                   else
                                                       if length(varargin)>3
                                                           if varargin{4}(1)=='q' && mod(str2double(varargin{4}(2)),1)==0
                                                               if str2double(varargin{4}(2))<=7 && str2double(varargin{4}(2))>=1 && length(varargin{4})<3
                                                                   rab2=str2double(varargin{4}(2)); 
                                                                   if length(varargin)>4
                                                                       if strcmp(varargin{5},'abs')
                                                                           abst=1;
                                                                           if length(varargin)>5
                                                                               disp('input after the end string (abs) will be ignored');
                                                                           end
                                                                       else
                                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                                       end 
                                                                   end
                                                               else
                                                                   error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                               end
                                                           else
                                                               if length(varargin)>3
                                                                   if strcmp(varargin{4},'abs')
                                                                       abst=1;
                                                                       if length(varargin)>4
                                                                           disp('input after the end string (abs) will be ignored');
                                                                       end
                                                                   else
                                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                                   end
                                                               end
                                                           end
                                                       end

                                                   end
                                               else
                                                   error('invalid input or improper order for the optional input arguments');
                                               end
                                           end
                                       else
                                           error('the given significance level must be a float number in (0,1)');
                                       end
                                   else 
                                       if length(varargin)>2
                                           if ischar(varargin{3})
                                               if varargin{3}(1)=='c'
                                                   if strcmp(varargin{3},'c1') || strcmp(varargin{3},'c2')
                                                       ctrig=str2double(varargin{3}(2));
                                                       if length(varargin)>3
                                                           if varargin{4}(1)=='q' && mod(str2double(varargin{4}(2)),1)==0
                                                               if str2double(varargin{4}(2))<=7 && str2double(varargin{4}(2))>=1 && length(varargin{4})<3
                                                                   rab2=str2double(varargin{4}(2)); 
                                                                   if length(varargin)>4
                                                                       if strcmp(varargin{5},'abs')
                                                                           abst=1;
                                                                           if length(varargin)>5
                                                                               disp('input after the end string (abs) will be ignored');
                                                                           end
                                                                       else
                                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                                       end
                                                                   end
                                                               else
                                                                   error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                               end
                                                           else
                                                               if length(varargin)>3
                                                                   if strcmp(varargin{4},'abs')
                                                                       abst=1;
                                                                       if length(varargin)>4
                                                                           disp('input after the end string (abs) will be ignored');
                                                                       end
                                                                   else
                                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                                   end
                                                               end
                                                           end
                                                       end
                                                   else
                                                       error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                                   end
                                               else
                                                   if length(varargin)>2
                                                       if varargin{3}(1)=='q' && mod(str2double(varargin{3}(2)),1)==0
                                                           if str2double(varargin{3}(2))<=7 && str2double(varargin{3}(2))>=1 && length(varargin{3})<3
                                                               rab2=str2double(varargin{3}(2)); 
                                                               if length(varargin)>3
                                                                   if strcmp(varargin{4},'abs')
                                                                       abst=1;
                                                                       if length(varargin)>4
                                                                           disp('input after the end string (abs) will be ignored');
                                                                       end
                                                                   else
                                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                                   end 
                                                               end
                                                           else
                                                               error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                           end
                                                       else
                                                           if length(varargin)>2
                                                               if strcmp(varargin{3},'abs')
                                                                   abst=1;
                                                                   if length(varargin)>3
                                                                       disp('input after the end string (abs) will be ignored');
                                                                   end
                                                               else
                                                                   error('unidentified input characters or improper order for the optional input arguments');
                                                               end 
                                                           end
                                                       end
                                                   end
                                               end
                                           else
                                               error('invalid input or improper order for the optional input arguments');
                                           end
                                       end
                                   end 
                               end 
                           else
                               if isfloat(varargin{2})
                                   error('min.perm must be a real integer number >=100');     
                               else
                                   error('invalid data type for the 2nd optional input argument'); 
                               end
                           end
                       end    
                   else
                       if isfloat(varargin{1}) 
                           if varargin{1}>0 && varargin{1}<1
                               thrvp=varargin{1};
                               if length(varargin)>1
                                   if ischar(varargin{2})
                                       if varargin{2}(1)=='c'
                                           if strcmp(varargin{2},'c1') || strcmp(varargin{2},'c2')
                                               ctrig=str2double(varargin{2}(2));
                                               if length(varargin)>2
                                                   if varargin{3}(1)=='q' && mod(str2double(varargin{3}(2)),1)==0
                                                       if str2double(varargin{3}(2))<=7 && str2double(varargin{3}(2))>=1 && length(varargin{3})<3
                                                           rab2=str2double(varargin{3}(2)); 
                                                           if length(varargin)>3
                                                               if strcmp(varargin{4},'abs') 
                                                                   abst=1;
                                                                   if length(varargin)>4
                                                                       disp('input after the end string (abs) will be ignored');
                                                                   end
                                                               else
                                                                   error('unidentified input characters or improper order for the optional input arguments');
                                                               end
                                                           end
                                                       else
                                                           error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                       end
                                                   else
                                                       if length(varargin)>2
                                                           if strcmp(varargin{3},'abs')
                                                               abst=1;
                                                               if length(varargin)>3
                                                                   disp('input after the end string (abs) will be ignored');
                                                               end
                                                           else
                                                               error('unidentified input characters or improper order for the optional input arguments');
                                                           end
                                                       end
                                                   end
                                               end
                                           else
                                               error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                           end
                                       else
                                           if length(varargin)>1
                                               if varargin{2}(1)=='q' && mod(str2double(varargin{2}(2)),1)==0
                                                   if str2double(varargin{2}(2))<=7 && str2double(varargin{2}(2))>=1 && length(varargin{2})<3
                                                       rab2=str2double(varargin{2}(2)); 
                                                       if length(varargin)>2
                                                           if strcmp(varargin{3},'abs')
                                                               abst=1;
                                                               if length(varargin)>3
                                                                   disp('input after the end string (abs) will be ignored');
                                                               end
                                                           else
                                                               error('unidentified input characters or improper order for the optional input arguments');
                                                           end 
                                                       end
                                                   else
                                                       error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                   end
                                               else 
                                                   if length(varargin)>1
                                                       if strcmp(varargin{2},'abs')
                                                           abst=1; 
                                                           if length(varargin)>2
                                                               disp('input after the end string (abs) will be ignored');
                                                           end
                                                       else
                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                       end 
                                                   end
                                               end
                                           end
                                       end
                                   else
                                       error('invalid input or improper order for the optional input arguments');
                                   end
                               end 
                           else 
                               error('the given significance level must be a real float number in (0,1)'); 
                           end
                       else 
                           if ischar(varargin{1})
                               if varargin{1}(1)=='c'
                                   if strcmp(varargin{1},'c1') || strcmp(varargin{1},'c2')
                                       ctrig=str2double(varargin{1}(2));
                                       if length(varargin)>1
                                           if varargin{2}(1)=='q' && mod(str2double(varargin{2}(2)),1)==0
                                               if str2double(varargin{2}(2))<=7 && str2double(varargin{2}(2))>=1 && length(varargin{2})<3
                                                   rab2=str2double(varargin{2}(2)); 
                                                   if length(varargin)>2
                                                       if strcmp(varargin{3},'abs')
                                                           abst=1;
                                                           if length(varargin)>3
                                                               disp('input after the end string (abs) will be ignored');
                                                           end
                                                       else
                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                       end
                                                   end
                                               else
                                                   error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                               end
                                           else
                                               if length(varargin)>1
                                                   if strcmp(varargin{2},'abs')
                                                       abst=1;
                                                       if length(varargin)>2
                                                           disp('input after the end string (abs) will be ignored');
                                                       end
                                                   else
                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                   end
                                               end
                                           end
                                       end
                                   else
                                       error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                   end
                               else
                                   if varargin{1}(1)=='q' && mod(str2double(varargin{1}(2)),1)==0
                                       if str2double(varargin{1}(2))<=7 && str2double(varargin{1}(2))>=1 && length(varargin{1})<3
                                           rab2=str2double(varargin{1}(2)); 
                                           if length(varargin)>1
                                               if strcmp(varargin{2},'abs')
                                                   abst=1;
                                                   if length(varargin)>2
                                                       disp('input after the end string (abs) will be ignored');
                                                   end
                                               else
                                                   error('unidentified input characters or improper order for the optional input arguments');
                                               end 
                                           end
                                       else
                                           error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                       end
                                   else
                                       if strcmp(varargin{1},'abs')
                                           abst=1;
                                           if length(varargin)>1
                                               disp('input after the end string (abs) will be ignored');
                                           end
                                       else
                                           error('unidentified input characters or improper order for the optional input arguments');
                                       end   
                                   end
                               end
                           else
                               error('invalid input or improper order for the optional input arguments');
                           end
                       end 
                   end         
            end
    end
    
    if length(varargin)>1 
        switch varargin{2}
            case 'Pearson'
               crtag=1;
            case 'Kendall'
               crtag=2;
            case 'Spearman'
               crtag=3;
            otherwise 
               if isfloat(varargin{2}) && mod(varargin{2},1)==0 && varargin{2}>=100 
                   if length(varargin)>2
                       if isfloat(varargin{3}) && mod(varargin{3},1)==0 && varargin{3}>=100
                           rnd=max(varargin{2},varargin{3});
                           rnde=min(varargin{2},varargin{3}); 
                           if length(varargin)>3
                               if isfloat(varargin{4})
                                   if varargin{4}>0 && varargin{4}<1 
                                       thrvp=varargin{4};
                                       if length(varargin)>4    
                                           if ischar(varargin{5})
                                               if varargin{5}(1)=='c'
                                                   if strcmp(varargin{5},'c1') || strcmp(varargin{5},'c2')
                                                       ctrig=str2double(varargin{5}(2));
                                                       if length(varargin)>5
                                                           if varargin{6}(1)=='q' && mod(str2double(varargin{6}(2)),1)==0
                                                               if str2double(varargin{6}(2))<=7 && str2double(varargin{6}(2))>=1 && length(varargin{6})<3 
                                                                   rab2=str2double(varargin{6}(2)); 
                                                                   if length(varargin)>6
                                                                       if strcmp(varargin{7},'abs') 
                                                                           abst=1;
                                                                           if length(varargin)>7
                                                                               disp('input after the end string (abs) will be ignored');
                                                                           end
                                                                       else
                                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                                       end 
                                                                   end
                                                               else
                                                                   error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                               end
                                                           else
                                                               if length(varargin)>5
                                                                   if strcmp(varargin{6},'abs')
                                                                       abst=1;
                                                                       if length(varargin)>6
                                                                           disp('input after the end string (abs) will be ignored');
                                                                       end
                                                                   else
                                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                                   end
                                                               end
                                                           end
                                                       end        
                                                   else
                                                       error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                                   end
                                               else
                                                   if length(varargin)>4
                                                       if varargin{5}(1)=='q' && mod(str2double(varargin{5}(2)),1)==0
                                                           if str2double(varargin{5}(2))<=7 && str2double(varargin{5}(2))>=1 && length(varargin{5})<3
                                                               rab2=str2double(varargin{5}(2)); 
                                                               if length(varargin)>5
                                                                   if strcmp(varargin{6},'abs')
                                                                       abst=1;
                                                                       if length(varargin)>6
                                                                           disp('input after the end string (abs) will be ignored');
                                                                       end
                                                                   else
                                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                                   end
                                                               end
                                                           else
                                                               error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                           end
                                                       else 
                                                           if length(varargin)>4
                                                               if strcmp(varargin{5},'abs')
                                                                   abst=1;
                                                                   if length(varargin)>5
                                                                       disp('input after the end string (abs) will be ignored');
                                                                   end
                                                               else
                                                                   error('unidentified input characters or improper order for the optional input arguments');
                                                               end
                                                           end
                                                       end
                                                   end
                                               end
                                           else
                                               error('invalid input or improper order for the optional input arguments');
                                           end
                                       end
                                   else
                                       error('the given significance level must be a float number in (0,1)');
                                   end
                               else 
                                   if length(varargin)>3
                                       if ischar(varargin{4})
                                           if varargin{4}(1)=='c'
                                               if strcmp(varargin{4},'c1') || strcmp(varargin{4},'c2')
                                                   ctrig=str2double(varargin{4}(2));
                                                   if length(varargin)>4
                                                       if varargin{5}(1)=='q' && mod(str2double(varargin{5}(2)),1)==0
                                                           if str2double(varargin{5}(2))<=7 && str2double(varargin{5}(2))>=1 && length(varargin{5})<3
                                                               rab2=str2double(varargin{5}(2)); 
                                                               if length(varargin)>5
                                                                   if strcmp(varargin{6},'abs')
                                                                       abst=1;
                                                                       if length(varargin)>6
                                                                           disp('input after the end string (abs) will be ignored');
                                                                       end
                                                                   else
                                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                                   end
                                                               end
                                                           else
                                                               error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                           end
                                                       else
                                                           if length(varargin)>4
                                                               if strcmp(varargin{5},'abs')
                                                                   abst=1;
                                                                   if length(varargin)>5
                                                                       disp('input after the end string (abs) will be ignored');
                                                                   end
                                                               else
                                                                   error('unidentified input characters or improper order for the optional input arguments');
                                                               end
                                                           end
                                                       end
                                                   end
                                               else
                                                   error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                               end
                                           else
                                               if length(varargin)>3
                                                   if varargin{4}(1)=='q' && mod(str2double(varargin{4}(2)),1)==0
                                                       if str2double(varargin{4}(2))<=7 && str2double(varargin{4}(2))>=1 && length(varargin{4})<3
                                                           rab2=str2double(varargin{4}(2)); 
                                                           if length(varargin)>4
                                                               if strcmp(varargin{5},'abs')
                                                                   abst=1;
                                                                   if length(varargin)>5
                                                                       disp('input after the end string (abs) will be ignored');
                                                                   end
                                                               else
                                                                   error('unidentified input characters or improper order for the optional input arguments');
                                                               end 
                                                           end
                                                       else
                                                           error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                       end
                                                   else
                                                       if length(varargin)>3
                                                           if strcmp(varargin{4},'abs')
                                                               abst=1;
                                                               if length(varargin)>4
                                                                   disp('input after the end string (abs) will be ignored');
                                                               end
                                                           else
                                                               error('unidentified input characters or improper order for the optional input arguments');
                                                           end 
                                                       end
                                                   end
                                               end
                                           end
                                       else
                                           error('invalid input or improper order for the optional input arguments');
                                       end
                                   end
                               end 
                           end 
                       else
                           if isfloat(varargin{1}) && mod(varargin{1},1)==0 && varargin{1}>=100 % nothing to do--keep rnde, rnd
                           else
                               if isfloat(varargin{3})
                                   error('max.perm must be a real integer number >=100');     
                               else
                                   error('invalid data type for the 2nd optional input argument'); 
                               end
                           end
                       end
                   end    
               else 
                   if isfloat(varargin{2}) 
                       if varargin{2}>0 && varargin{2}<1
                           thrvp=varargin{2};
                           if length(varargin)>2
                               if ischar(varargin{3})
                                   if varargin{3}(1)=='c'
                                       if strcmp(varargin{3},'c1') || strcmp(varargin{3},'c2')
                                           ctrig=str2double(varargin{3}(2));
                                           if length(varargin)>3
                                               if varargin{4}(1)=='q' && mod(str2double(varargin{4}(2)),1)==0
                                                   if str2double(varargin{4}(2))<=7 && str2double(varargin{4}(2))>=1 && length(varargin{4})<3
                                                       rab2=str2double(varargin{4}(2)); 
                                                       if length(varargin)>4
                                                           if strcmp(varargin{5},'abs') 
                                                               abst=1;
                                                               if length(varargin)>5
                                                                   disp('input after the end string (abs) will be ignored');
                                                               end
                                                           else
                                                               error('unidentified input characters or improper order for the optional input arguments');
                                                           end
                                                       end
                                                   else
                                                       error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                   end
                                               else
                                                   if length(varargin)>3
                                                       if strcmp(varargin{4},'abs')
                                                           abst=1;
                                                           if length(varargin)>4
                                                               disp('input after the end string (abs) will be ignored');
                                                           end
                                                       else
                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                       end
                                                   end
                                               end
                                           end
                                       else
                                           error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                       end
                                   else
                                       if length(varargin)>2
                                           if varargin{3}(1)=='q' && mod(str2double(varargin{3}(2)),1)==0
                                               if str2double(varargin{3}(2))<=7 && str2double(varargin{3}(2))>=1 && length(varargin{3})<3
                                                   rab2=str2double(varargin{3}(2)); 
                                                   if length(varargin)>3
                                                       if strcmp(varargin{4},'abs')
                                                           abst=1;
                                                           if length(varargin)>4
                                                               disp('input after the end string (abs) will be ignored');
                                                           end
                                                       else
                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                       end 
                                                   end
                                               else
                                                   error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                               end
                                           else
                                               if length(varargin)>2
                                                   if strcmp(varargin{3},'abs')
                                                       abst=1;
                                                       if length(varargin)>3
                                                           disp('input after the end string (abs) will be ignored');
                                                       end
                                                   else
                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                   end 
                                               end
                                           end
                                       end
                                   end
                               else
                                   error('invalid input or improper order for the optional input arguments');
                               end
                           end 
                       else 
                           error('the given significance level must be a real float number in (0,1)'); 
                       end
                   else 
                       if ischar(varargin{2})
                           if varargin{2}(1)=='c'
                               if strcmp(varargin{2},'c1') || strcmp(varargin{2},'c2')
                                   ctrig=str2double(varargin{2}(2));
                                   if length(varargin)>2
                                       if varargin{3}(1)=='q' && mod(str2double(varargin{3}(2)),1)==0
                                           if str2double(varargin{3}(2))<=7 && str2double(varargin{3}(2))>=1 && length(varargin{3})<3
                                               rab2=str2double(varargin{3}(2)); 
                                               if length(varargin)>3
                                                   if strcmp(varargin{4},'abs')
                                                       abst=1;
                                                       if length(varargin)>4
                                                           disp('input after the end string (abs) will be ignored');
                                                       end
                                                   else
                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                   end
                                               end
                                           else
                                               error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                           end
                                       else
                                           if length(varargin)>2
                                               if strcmp(varargin{3},'abs')
                                                   abst=1;
                                                   if length(varargin)>3
                                                       disp('input after the end string (abs) will be ignored');
                                                   end
                                               else
                                                   error('unidentified input characters or improper order for the optional input arguments');
                                               end
                                           end
                                       end
                                   end
                               else
                                   error('invalid input for individual cluster quality measure (must be c1 or c2)');
                               end
                           else
                               if varargin{2}(1)=='q' && mod(str2double(varargin{2}(2)),1)==0
                                   if str2double(varargin{2}(2))<=7 && str2double(varargin{2}(2))>=1 && length(varargin{2})<3
                                       rab2=str2double(varargin{2}(2)); 
                                       if length(varargin)>2
                                           if strcmp(varargin{3},'abs')
                                               abst=1;
                                               if length(varargin)>3
                                                   disp('input after the end string (abs) will be ignored');
                                               end
                                           else
                                               error('unidentified input characters or improper order for the optional input arguments');
                                           end
                                       end
                                   else
                                       error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                   end
                               else
                                   if strcmp(varargin{2},'abs')
                                       abst=1;
                                       if length(varargin)>2
                                           disp('input after the end string (abs) will be ignored');
                                       end
                                   else
                                       error('unidentified input characters or improper order for the optional input arguments');
                                   end   
                               end
                           end
                       else
                           error('invalid input or improper order for the optional input arguments');
                       end
                   end 
               end 
        end
    end
    
    if length(varargin)>2 
          if isfloat(varargin{3}) && mod(varargin{3},1)==0 && varargin{3}>=100 
               if length(varargin)>3
                   if isfloat(varargin{4}) && mod(varargin{4},1)==0 && varargin{4}>=100
                       rnd=max(varargin{3},varargin{4});
                       rnde=min(varargin{3},varargin{4}); 
                       if length(varargin)>4
                           if isfloat(varargin{5})
                               if varargin{5}>0 && varargin{5}<1 
                                   thrvp=varargin{5};
                                   if length(varargin)>5    
                                       if ischar(varargin{6})
                                           if varargin{6}(1)=='c'
                                               if strcmp(varargin{6},'c1') || strcmp(varargin{6},'c2')
                                                   ctrig=str2double(varargin{6}(2));
                                                   if length(varargin)>6
                                                       if varargin{7}(1)=='q' && mod(str2double(varargin{7}(2)),1)==0
                                                           if str2double(varargin{7}(2))<=7 && str2double(varargin{7}(2))>=1 && length(varargin{7})<3 
                                                               rab2=str2double(varargin{7}(2)); 
                                                               if length(varargin)>7
                                                                   if strcmp(varargin{8},'abs') 
                                                                       abst=1;
                                                                       if length(varargin)>8
                                                                           disp('input after the end string (abs) will be ignored');
                                                                       end
                                                                   else
                                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                                   end
                                                               end
                                                           else
                                                               error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                           end
                                                       else
                                                           if length(varargin)>6
                                                               if strcmp(varargin{7},'abs')
                                                                   abst=1;
                                                                   if length(varargin)>7
                                                                       disp('input after the end string (abs) will be ignored');
                                                                   end
                                                               else
                                                                   error('unidentified input characters or improper order for the optional input arguments');
                                                               end
                                                           end
                                                       end
                                                   end        
                                               else
                                                   error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                               end
                                           else
                                               if length(varargin)>5
                                                   if varargin{6}(1)=='q' && mod(str2double(varargin{6}(2)),1)==0
                                                       if str2double(varargin{6}(2))<=7 && str2double(varargin{6}(2))>=1 && length(varargin{6})<3
                                                           rab2=str2double(varargin{6}(2)); 
                                                           if length(varargin)>6
                                                               if strcmp(varargin{7},'abs')
                                                                   abst=1;
                                                                   if length(varargin)>7
                                                                       disp('input after the end string (abs) will be ignored');
                                                                   end
                                                               else
                                                                   error('unidentified input characters or improper order for the optional input arguments');
                                                               end
                                                           end
                                                       else
                                                           error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                       end
                                                   else
                                                       if length(varargin)>5
                                                           if strcmp(varargin{6},'abs')
                                                               abst=1;
                                                               if length(varargin)>6
                                                                   disp('input after the end string (abs) will be ignored');
                                                               end
                                                           else
                                                               error('unidentified input characters or improper order for the optional input arguments');
                                                           end
                                                       end
                                                   end
                                               end

                                           end
                                       else
                                           error('invalid input or improper order for the optional input arguments');
                                       end
                                   end
                               else
                                   error('the given significance level must be a float number in (0,1)');
                               end
                           else 
                               if length(varargin)>4
                                   if ischar(varargin{5})
                                       if varargin{5}(1)=='c'
                                           if strcmp(varargin{5},'c1') || strcmp(varargin{5},'c2')
                                               ctrig=str2double(varargin{5}(2));
                                               if length(varargin)>5
                                                   if varargin{6}(1)=='q' && mod(str2double(varargin{6}(2)),1)==0
                                                       if str2double(varargin{6}(2))<=7 && str2double(varargin{6}(2))>=1 && length(varargin{6})<3
                                                           rab2=str2double(varargin{6}(2)); 
                                                           if length(varargin)>6
                                                               if strcmp(varargin{7},'abs')
                                                                   abst=1;
                                                                   if length(varargin)>7
                                                                       disp('input after the end string (abs) will be ignored');
                                                                   end
                                                               else
                                                                   error('unidentified input characters or improper order for the optional input arguments');
                                                               end
                                                           end
                                                       else
                                                           error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                       end
                                                   else
                                                       if length(varargin)>5
                                                           if strcmp(varargin{6},'abs')
                                                               abst=1;
                                                               if length(varargin)>6
                                                                   disp('input after the end string (abs) will be ignored');
                                                               end
                                                           else
                                                               error('unidentified input characters or improper order for the optional input arguments');
                                                           end
                                                       end
                                                   end
                                               end
                                           else
                                               error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                           end
                                       else
                                           if length(varargin)>4
                                               if varargin{5}(1)=='q' && mod(str2double(varargin{5}(2)),1)==0
                                                   if str2double(varargin{5}(2))<=7 && str2double(varargin{5}(2))>=1 && length(varargin{5})<3
                                                       rab2=str2double(varargin{5}(2)); 
                                                       if length(varargin)>5
                                                           if strcmp(varargin{6},'abs')
                                                               abst=1;
                                                               if length(varargin)>6
                                                                   disp('input after the end string (abs) will be ignored');
                                                               end
                                                           else
                                                               error('unidentified input characters or improper order for the optional input arguments');
                                                           end
                                                       end
                                                   else
                                                       error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                                   end
                                               else
                                                   if length(varargin)>4
                                                       if strcmp(varargin{5},'abs')
                                                           abst=1;
                                                           if length(varargin)>5
                                                               disp('input after the end string (abs) will be ignored');
                                                           end
                                                       else
                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                       end 
                                                   end
                                               end
                                           end
                                       end
                                   else
                                       error('invalid input or improper orde for the optional input arguments');
                                   end
                               end
                           end
                       end 
                   else
                       if isfloat(varargin{2}) && mod(varargin{2},1)==0 && varargin{2}>=100
                       else
                           if isfloat(varargin{4})
                               error('min.perm must be a real integer number >=100');     
                           else
                               error('invalid data type for the 2nd optional input argument'); 
                           end
                       end
                   end
               end 
          else
               if isfloat(varargin{3}) 
                   if varargin{3}>0 && varargin{3}<1
                       thrvp=varargin{3};
                       if length(varargin)>3
                           if ischar(varargin{4})
                               if varargin{4}(1)=='c'
                                   if strcmp(varargin{4},'c1') || strcmp(varargin{4},'c2')
                                       ctrig=str2double(varargin{4}(2));
                                       if length(varargin)>4
                                           if varargin{5}(1)=='q' && mod(str2double(varargin{5}(2)),1)==0
                                               if str2double(varargin{5}(2))<=7 && str2double(varargin{5}(2))>=1 && length(varargin{5})<3
                                                   rab2=str2double(varargin{5}(2)); 
                                                   if length(varargin)>5
                                                       if strcmp(varargin{6},'abs') 
                                                           abst=1;
                                                           if length(varargin)>6
                                                               disp('input after the end string (abs) will be ignored');
                                                           end
                                                       else
                                                           error('unidentified input characters or improper order for the optional input arguments');
                                                       end
                                                   end
                                               else
                                                   error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                               end
                                           else
                                               if length(varargin)>4
                                                   if strcmp(varargin{5},'abs')
                                                       abst=1;
                                                       if length(varargin)>5
                                                           disp('input after the end string (abs) will be ignored');
                                                       end
                                                   else
                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                   end
                                               end
                                           end
                                       end
                                   else
                                       error('invalid input for individual cluster quality measure (must be c1 or c2)');
                                   end
                               else
                                   if length(varargin)>3
                                       if varargin{4}(1)=='q' && mod(str2double(varargin{4}(2)),1)==0
                                           if str2double(varargin{4}(2))<=7 && str2double(varargin{4}(2))>=1 && length(varargin{4})<3
                                               rab2=str2double(varargin{4}(2)); 
                                               if length(varargin)>4
                                                   if strcmp(varargin{5},'abs')
                                                       abst=1;
                                                       if length(varargin)>5
                                                           disp('input after the end string (abs) will be ignored');
                                                       end
                                                   else
                                                       error('unidentified input characters or improper order for the optional input arguments');
                                                   end
                                               end
                                           else
                                               error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                           end
                                       else 
                                           if length(varargin)>3
                                               if strcmp(varargin{4},'abs')
                                                   abst=1;
                                                   if length(varargin)>4
                                                       disp('input after the end string (abs) will be ignored');
                                                   end
                                               else
                                                   error('unidentified input characters or improper order for the optional input arguments');
                                               end 
                                           end
                                       end
                                   end
                               end
                           else
                               error('invalid input or improper orde for the optional input arguments');
                           end
                       end
                   else 
                       error('the given significance level must be a real float number in (0,1)'); 
                   end
               else 
                   if ischar(varargin{3})
                       if varargin{3}(1)=='c'
                           if strcmp(varargin{3},'c1') || strcmp(varargin{3},'c2')
                               ctrig=str2double(varargin{3}(2));
                               if length(varargin)>2
                                   if varargin{4}(1)=='q' && mod(str2double(varargin{4}(2)),1)==0
                                       if str2double(varargin{4}(2))<=7 && str2double(varargin{4}(2))>=1 && length(varargin{4})<3
                                           rab2=str2double(varargin{4}(2)); 
                                           if length(varargin)>4
                                               if strcmp(varargin{5},'abs')
                                                   abst=1;
                                                   if length(varargin)>5
                                                       disp('input after the end string (abs) will be ignored');
                                                   end
                                               else
                                                   error('unidentified input characters or improper order for the optional input arguments');
                                               end
                                           end
                                       else
                                           error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                                       end
                                   else
                                       if length(varargin)>3
                                           if strcmp(varargin{4},'abs')
                                               abst=1;
                                               if length(varargin)>4
                                                   disp('input after the end string (abs) will be ignored');
                                               end
                                           else
                                               error('unidentified input characters or improper order for the optional input arguments');
                                           end
                                       end
                                   end
                               end
                           else
                               error('invalid input for individual cluster quality measure (must be c1 or c2)');
                           end
                       else
                           if varargin{3}(1)=='q' && mod(str2double(varargin{3}(2)),1)==0
                               if str2double(varargin{3}(2))<=7 && str2double(varargin{3}(2))>=1 && length(varargin{3})<3
                                   rab2=str2double(varargin{3}(2)); 
                                   if length(varargin)>3
                                       if strcmp(varargin{4},'abs')
                                           abst=1;
                                           if length(varargin)>4
                                               disp('input after the end string (abs) will be ignored');
                                           end
                                       else
                                           error('unidentified input characters or improper order for the optional input arguments');
                                       end
                                   end
                               else
                                   error('invalid input for overall cluster quality measure (must be q1, q2, ..., q7)');
                               end
                           else
                               if strcmp(varargin{3},'abs')
                                   abst=1;
                                   if length(varargin)>3
                                       disp('input after the end string (abs) will be ignored');
                                   end
                               else
                                   error('unidentified input characters or improper order for the optional input arguments');
                               end   
                           end
                       end
                   else
                       error('invalid input or improper orde for the optional input arguments');
                   end
               end 
          end
    end
 end 
                                        
end     
%-------------------------------------------------------------------------%
% End of nested function: checkinputs
%-------------------------------------------------------------------------%

end % end of C3D