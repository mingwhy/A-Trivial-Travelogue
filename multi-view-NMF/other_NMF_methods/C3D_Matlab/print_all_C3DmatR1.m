function [ds,pv,itr_flag,rab,pvrab,itr_flagrab,consout]=validate_trans(data_type,D,V,coln,td1,td2,rnd,rnde,thrvp,ctrig,rab2,abst)
% General process to extract the conditions (data matrices) where the cluster(s) are
% present and to compute their corresponding significance (by computing a
% p-value)
% Input (Mandatory)
%   data_type: co-expression or expression
%   D: multiple real input data sets, the input data sets must have an identical 
%   ROW size 
%   V: Common factor
%   coln: candidate vector
%   td1, td2: nodes positions of the cluster 
% Input (Optional)
%    rnd      = max.perm
%    rnde     = min.perm
%    thrvp    = p.threshold
%    ctrig    = individual.cluster.quality
%    rab2     = overall.cluster.quality
%    abst     = norm
%
% Output
%   ds: inidvidual cluster quality 
%   pv: individual p-value 
%   itr_flag: permutations (for pv)
%   rab: overall cluster quality
%   pvrab: overall p-value
%   itr_flagrab:  permutations (for pvrab)
%   consout: conditions

% Written by Xiaolin Xiao, 2012


if nargin<11
    error(message('TooFewInputs'));
end
if (nargout > 7)
    error(message('MATLAB:clustervalidation:TooManyOutputs'));
end

% essential input arguments: routine check (data type)
if isreal(V)==0
    error('common factor has complex numbers');
end
if mod(td1,1) || mod(td2,1)
    error('node number must be an integer');
end

% essential input arguments: other check
if td1>=td2
    error('can not locate the target cluster');
else
    intv=td2-td1+1;
end

cond=length(D);
if cond<2
    error('number of conditions is less than 2');
end
all=1:1:cond;

drsize=nan(1,cond);
real_tag=nan(1,cond);
for i=1:cond
    drsize(i)=size(D{i},1);
    real_tag(i)=isreal(D{i});
end
if isequal(drsize, drsize(1).*ones(1,cond))==0
    error(message('MATLAB:hogsvd:MatrixRowMismatch'))
end
if isequal(real_tag,ones(1,cond))==0
    error('Complex numbers exist');
end

if isreal(V)==0
    error('Incorrect Common Factor');
end

if isequal(size(D{1},1),length(V))==0
    error('Common factor and input data matrices mismatch');
end

if coln>length(V) || mod(coln,1)~=0 
    error('Incorrect input of the eigen-vector');
end

N=length(V);
if td1>N || td2>N
    error('cluster not found');
end

% post-argument checking
one=1; 
sv=1e-16; 


rsize=nan(cond,1);
csize=nan(cond,1);
if data_type~=1 % nonsymetric data matrices
    Dc=cell(cond,1);
    for i=1:cond
        Dc{i}=D{i}*D{i}'; 
        rsize(i)=size(Dc{i},1);
        csize(i)=size(Dc{i},2);
        if rsize(i)~=mean(rsize(1:i))
            error('input data set %d row size mismatch',i); 
        end
        if csize(i)~=mean(csize(1:i))
            error('input data set %d column size mismatch',i);
        end
        % not tolerate Nan(missing number), Inf and complex numbers
        if sum(sum(isnan(Dc{i})))>0
           error('Nan in data (symmetric data matrix generated), re-try with different arguments');
        end
        if sum(sum(isinf(Dc{i})))>0
            error('Inf in data (symmetric data matrix generated), re-try with different arguments');
        end
        if ~isreal(Dc{i})
            error('Complex number existin data (symmetric data matrix generated), re-try with different arguments');
        end
    end
else % co-expression
    Dc=D; 
    for i=1:cond
        rsize(i)=size(Dc{i},1);
        csize(i)=size(Dc{i},2);
        if rsize(i)~=mean(rsize(1:i))
            error('input data set %d row size mismatch',i);
        end
        if csize(i)~=mean(csize(1:i))
            error('input data set %d column size mismatch',i);
        end
    end     
end
clear D;
% same row size
if isequal(length(V),length(Dc{1}))==0
    error('Size mismatch: Common factor and input data matrices');
end
   
if abst==1
    Dc=cellfun(@abs,Dc,'un',0);
end
    
[~,b] = sort(V(:,coln));
Ds=cell(cond,1);
for i=1:cond
    Ds{i}=Dc{i}(b,b);
end


if ctrig==2 
    di=nan(1,cond);
    do=nan(1,cond);
    ds=nan(1,cond);
    for i=1:cond
        di(i)=sum(sum(Ds{i}(td1:td2,td1:td2),2));
        do(i)=(sum(sum(Ds{i},2))-di(i))./(2*(N-intv));
    end    
    di=di./(2*intv);
    do(do==0)=sv;
    ds=di./do; 
else
    if ctrig==1
        ds=nan(1,cond);
        for i=1:cond
            ds(i)=mean(mean(Ds{i}(td1:td2,td1:td2)));
        end   
    else
        error('Options for individual cluster quality must be 1 or 2!'); 
    end
end

if sum(isnan(ds))>0 || sum(isinf(ds))>0
    error('Mistake: individual cluster quality computation');
end

% sub-step 1--individual permutation test, min.perm (rnde)
disp('Validation sub-step 1--permuting each input data set individually...');

dse=nan(cond,rnde);
for i=1:rnde
    pr=randperm(N); 
    for j=1:cond
        Ds{j}=Dc{j}(pr,pr);
    end 
    if ctrig==2 
        for j=1:cond
            di(j)=sum(sum(Ds{j}(1:intv,1:intv),2));
            do(j)=(sum(sum(Ds{j},2))-di(j))./(2*(N-intv));
        end
        di=di./(2*intv);
        do(do==0)=sv;
        dse(:,i)=di./do;
    else
        for j=1:cond
            dse(j,i)=mean(mean(Ds{j}(1:intv,1:intv)));
        end
    end   
end

if sum(sum(isnan(dse)))>0 || sum(sum(isinf(dse)))>0
    error('error in computing individual cluster quality samples');
end

pv=nan(1,cond);
for i=1:cond
    pv(i)=length(find(dse(i,:)>=ds(i)))/rnde; 
    if isnan(pv(i))==1 || isinf(pv(i))==1
        error('Validation sub-step 1 (min.perm): error in computing this p-value for the data set');
    end
end

itr_flag=nan(1,cond);
itr_flag(:)=rnde;

% sub-step 1---individual incremental permutation test (to max.perm rnd)
if rnd>rnde
    
    thr=thrvp; 
    cons=find(pv<=thr);
    cons=sort(cons);
    conslen=length(cons);
    
    if conslen >0 
        disp('incremental permutation begins...');
        dss=nan(conslen,rnd);
        dss(:,1:rnde)=dse(cons,:);
        dset=nan(cond,rnd);
        dset(:,1:rnde)=dse; 
        clear dse;      
        i=rnde+1;
        while i<=rnd
            pr=randperm(N);
            for j=1:conslen
                Ds{cons(j)}=Dc{cons(j)}(pr,pr);
            end    
            di=nan(1,conslen);
            do=nan(1,conslen);   
            if ctrig==2 
                for j=1:conslen
                    di(j)=sum(sum(Ds{cons(j)}(1:intv,1:intv),2));
                    do(j)=(sum(sum(Ds{cons(j)},2))-di(j))./(2*(N-intv));
                end              
                di=di./(2*intv);
                do(do==0)=sv;
                dss(:,i)=di./do;
            else
                for j=1:conslen
                    dss(j,i)=mean(mean(Ds{cons(j)}(1:intv,1:intv)));
                end
            end
            dset(cons,i)=dss(:,i);
            for j=1:conslen
                itr_flag(cons(j))=i;
                pv(cons(j))=length(find(dss(j,1:i)>=ds(cons(j))))/i;
                if isnan(pv(cons(j)))==1 || isinf(pv(cons(j)))==1
                    error('Validation sub-step 1 incremental permutation: error in computing the p-value for data set %d\n',cons(j));
                end    
            end
          
            [~,ia,~]=intersect(cons,cons(pv(cons)>thr));
            dss(ia,:)=[];
            cons=setdiff(cons,cons(pv(cons)>thr)); cons=sort(cons); 
           
            if conslen > length(cons)
                conslen=length(cons);
            else
                if conslen < length(cons)
                    error('G+ size mismatch, return');
                end
            end
            
            if conslen>0 
                i=i+1;
            else 
                i=rnd+1; 
            end       
        end      
    end      
end
    
clear dss di do;  


disp('Validation Step 1 Done');
disp('Here are the output from Validation Step 1:');
for i=1:cond
    if isnan(pv(i))==1 || isnan(itr_flag(i))==1
        error('Validation sub-step 1 incremental permutation: error in computing itr_flag for data set %d\n',i);
    else
        if isinf(pv(i))==0 && isinf(itr_flag(i))==0
            fprintf('Data set %d: the individual p-value is %f by randomizing the data %d times\n',i,pv(i),itr_flag(i));
        else
            error('Validation sub-step 1 incremental permutation: error in computing pv or itr_flag');
        end
    end
end

cons=sort(cons); 
disp('The cluster presents in data sets (G+):'); disp(cons);


% sub-step 2 -- incremental permutation test
disp('Validation Step 2---validating the significance of the cluster in data sets G+');

if length(find(pv<=thr)) ~= conslen || length(cons) ~= length(find(pv<=thr)) 
    error('mistake in Validation Step 1, check...'); 
else    
    if conslen>0
        if isequal(length(cons),conslen)==0
            error('length of G+ mismatch');
        end
        consout=cons; 
        constop=setdiff(all,cons);
        constop=sort(constop); 
        stoplen=length(constop); 
        
        if stoplen==0 && conslen==1 
            error('mistake in Validation sub-step 1, number of input data sets mismatch!');
        end
        
        if stoplen<0 
            error('mistake in Validation sub-step 1, number of input data sets mismatch!');
        end
        
        if (conslen+stoplen)~=cond
            error('size of G mismatch');
        end
         
        if sum(isnan(ds))>0 || sum(isinf(ds))>0
            error('Mistake: individual cluster quality computation');
        end

        
        switch rab2 
           case 1 
               if isempty(constop) 
                   det=one;
               else 
                   det=ds(constop);
                   det(det==0)=sv; 
               end
               rab=prod(ds(cons))./(prod(det));  
               clear det;
           case 2  
               if isempty(constop)
                   det=one;
               else
                   det=sum(ds(constop));
                   det(det==0)=sv;
               end
               rab=sum(ds(cons))./det;
               clear det;
           case 3  
               if isempty(constop) 
                   det=one;
               else
                   det=mean(ds(constop));
                   det(det==0)=sv; 
               end
               rab=mean(ds(cons))./det;
               clear det;
           case 4  
               if isempty(constop) 
                   det=one;
               else
                   det=geomean(ds(constop));
                   det(det==0)=sv;
               end
               rab=geomean(ds(cons))./det;
               clear det;
           case 5 
                if conslen>1
                    rab1s=conslen*(conslen-1)*prod(ds(cons));
                else
                    if conslen==1
                        rab1s=prod(ds(cons)); 
                    else 
                        error('data sets G+ size mismatch');
                    end
                end
                det=ds(constop);
                det(det==0)=sv;
                if isempty(constop) 
                    rab2s=one; 
                else   
                    if stoplen>1
                        rab2s=stoplen*(stoplen-1)*prod(det);          
                    else
                        rab2s=prod(det);
                    end 
                end
                rab=rab1s./rab2s;  
                clear rab1s rab2s det;
           case 6    
               det=ds;
               det(det==0)=sv; 

                if conslen>1
                    rab1s=0;
                    for i=1:conslen-1
                       for j=i+1:conslen
                           rab1s=ds(cons(i))./det(cons(j))+ds(cons(j))./det(cons(i))+rab1s; 
                       end
                    end
                    rab1s=conslen*(conslen-1)*rab1s; 
                else 
                     if conslen==1
                         rab1s=ds(cons); 
                     else
                         error('data sets G+ size mismatch');
                     end
                end
                if isempty(constop) 
                    rab2s=one;
                else
                    if stoplen>1
                        rab2s=0;
                        for i=1:stoplen-1
                           for j=i+1:stoplen
                               rab2s=ds(constop(i))./det(constop(j))+ds(constop(j))./det(constop(i))+rab2s;
                           end
                        end
                        rab2s(rab2s==0)=sv;
                        rab2s=stoplen*(stoplen-1)*rab2s; 
                    else
                         if stoplen==1
                             rab2s=det(constop);
                         else
                             error('data sets G- size mismatch');
                         end
                    end
                end
                rab=rab1s./rab2s; 
                clear rab1s rab2s det;
           case 7 
               dets=ds;
               dets(dets==0)=sv; 
               if conslen>1 
                   rab1s=0;
                   for i=1:conslen-1
                       for j=i+1:conslen
                           rab1s=ds(cons(i))./dets(cons(j))+ds(cons(j))./dets(cons(i))+rab1s;
                       end
                   end 
                   rab1s(rab1s==0)=sv;
               else 
                    if conslen==1
                         rab1s=dets(cons);   
                    else 
                         error('data sets G+ size mismatch');
                    end
               end
               if isempty(constop)
                   det=one;
                   rab2s=one;
               else 
                   if stoplen>1
                       rab2s=0;
                       for i=1:stoplen-1
                           for j=i+1:stoplen
                               rab2s=ds(constop(i))./dets(constop(j))+ds(constop(j))./dets(constop(i))+rab2s;
                           end
                       end
                       det=sum(ds(constop)); 
                       det(det==0)=sv;
                   else 
                       if stoplen==1
                           rab2s=ds(constop);
                           det=dets(constop); 
                       else 
                           error('data sets G- size mismatch');
                       end
                   end
                   
               end
               rab1s=det.*rab1s;              
               rab2s=sum(ds(cons)).*rab2s; 
               rab=rab2s./rab1s;
               clear rab1s rab2s det dets;
            otherwise
               error('Invalid input (1-7) for cluster quality');
       end        
        
        if isnan(rab)==0 && isinf(rab)==0
            fprintf('the original overall cluster quality is %f\n',rab); 
        else
            error('Mistake: overall cluster quality computation');
        end
        
        % sub-step 2, min.perm
        if ismember('dset',who)==1
            dse=dset(:,1:rnde); 
        else
            if ismember('dse',who)==0
                error('mistake in Validation sub-step 1 cluster validation'); 
            end
        end
        
       if sum(sum(isnan(dse(cons,:))))>0 || sum(sum(isinf(dse(cons,:))))>0
            error('Mistake: Validation sub-step 1 (min.perm)');
       end

       if size(dse,2)~=rnde
           error('min.perm number mismatch');
       end

       rabse=nan(1,rnde);      
       switch rab2
           case 1
               if isempty(constop)
                   det=one.*ones(1,rnde);
               else
                   det=dse(constop,:);  
                   det(det==0)=sv; 
                   if stoplen>1
                       det=prod(det);
                   end
               end
               if conslen>1
                   rabse=prod(dse(cons,:))./det;
               else
                   if conslen==1
                       rabse=dse(cons,:)./det;
                   else 
                       error('data sets G+ size mismatch');
                   end
               end
               clear det;
           case 2 
               if isempty(constop)
                   det=one.*ones(1,rnde);
               else
                   if stoplen==1
                       det=dse(constop,:);
                   else
                       det=sum(dse(constop,:));
                   end
                   det(det==0)=sv; 
               end
               if conslen>1
                   rabse=sum(dse(cons,:))./det;
               else
                   if conslen==1
                       rabse=dse(cons,:)./det;
                   else
                       error('data sets G+ size mismatch');
                   end
               end
               clear det;
           case 3
               if isempty(constop)
                   det=one.*ones(1,rnde);
               else
                   if stoplen==1
                       det=dse(constop,:);
                   else
                       det=mean(dse(constop,:));
                   end
                   det(det==0)=sv; 
               end
               if conslen>1
                   rabse=mean(dse(cons,:))./det;
               else
                   if conslen==1
                       rabse=dse(cons,:)./det;
                   else
                        error('data sets G+ size mismatch');
                   end
               end  
               clear det;
           case 4
               if isempty(constop) 
                   det=one.*ones(1,rnde);
               else
                   if stoplen==1
                       det=dse(constop,:);
                   else
                       det=geomean(dse(constop,:));
                   end
                   det(det==0)=sv; 
               end
               if conslen>1
                   rabse=geomean(dse(cons,:))./det;
               else
                   if conslen==1
                       rabse=dse(cons,:)./det;
                   else
                        error('data sets G+ size mismatch');
                   end
               end  
               clear det;                   
           case 5
               rabse1=nan(1,rnde);
               rabse2=nan(1,rnde);
               if conslen>1
                   rabse1=conslen*(conslen-1)*prod(dse(cons,:));
               else
                   if conslen==1
                       rabse1=dse(cons,:);
                   else
                       error('data sets G+ size mismatch');
                   end
               end
               det=dse(constop,:);
               det(det==0)=sv; 
               if isempty(constop)
                    rabse2=one.*ones(1,rnde); 
               else
                   if stoplen>1 
                       rabse2=stoplen*(stoplen-1)*prod(det);
                   else
                       if stoplen==1
                           rabse2=det;
                       end
                   end
               end
               rabse=rabse1./rabse2;
               clear rabse1 rabse2 det; 
           case 6
               det=dse;
               det(det==0)=sv;
               rabse1=nan(1,rnde);
               rabse2=nan(1,rnde);
               if conslen>1
                   rabse1=zeros(1,rnde);
                   for k=1:conslen-1
                       for j=k+1:conslen 
                           rabse1=dse(cons(k),:)./det(cons(j),:)+dse(cons(j),:)./det(cons(k),:)+rabse1;
                       end
                   end
                   rabse1=conslen*(conslen-1)*rabse1;
               else
                   if conslen==1
                       rabse1=dse(cons,:);
                   else
                       error('data sets G+ size mismatch');
                   end      
               end
               if isempty(constop)
                   rabse2=one.*ones(1,rnde);
               else    
                   if stoplen>1
                       rabse2=zeros(1,rnde);
                       for k=1:stoplen-1
                           for j=k+1:stoplen
                               rabse2=dse(constop(k),:)./det(constop(j),:)+dse(constop(j),:)./det(constop(k),:)+rabse2;
                           end
                       end
                       rabse2(rabse2==0)=sv;
                       rabse2=stoplen*(stoplen-1)*rabse2;
                   else
                       if stoplen==1
                            rabse2=det(constop,:); 
                       else 
                            error('data sets G- size mismatch');
                       end             
                   end
               end
               rabse=rabse1./rabse2; 
               clear rabse1 rabse2 det;           
           case 7
               det=dse;
               det(det==0)=sv;
               rabse1=nan(1,rnde);
               rabse2=nan(1,rnde);
               if conslen>1
                   det1=sum(dse(cons,:));
                   rabse1=zeros(1,rnde);
                   for k=1:conslen-1
                       for j=k+1:conslen
                           rabse1=dse(cons(k),:)./det(cons(j),:)+dse(cons(j),:)./det(cons(k),:)+rabse1;
                       end
                   end
                   rabse1(rabse1==0)=sv;
               else
                   if conslen==1
                       det1=dse(cons,:);
                       rabse1=det(cons,:); 
                   else
                        error('data sets G+ size mismatch');
                   end        
               end
               if isempty(constop) 
                   det2=one.*ones(1,rnde); 
                   rabse2=one.*ones(1,rnde); 
               else       
                   if stoplen>1
                       det2=sum(dse(constop,:)); 
                       det2(det2==0)=sv;
                       rabse2=zeros(1,rnde);
                       for k=1:stoplen-1
                           for j=k+1:stoplen
                               rabse2=dse(constop(k),:)./det(constop(j),:)+dse(constop(j),:)./det(constop(k),:)+rabse2;
                           end
                       end
                       
                   else
                       if stoplen==1
                           det2=det(constop,:); 
                           rabse2=dse(constop,:);
                       else 
                           error('data sets G- size mismatch');
                       end
                   end     
               end
               rabse2=det1.*rabse2;           
               rabse1=det2.*rabse1; 
               rabse=rabse2./rabse1;
               clear rabse1 rabse2 det1 det2 det;   
           otherwise
               error('Invalid input (1-7) for cluster quality');               
       end     
        
        
        if sum(isnan(rabse))>0 || sum(isinf(rabse))>0
            error('min.perm test mistake');
        end
        
        pvrab=length(find(rabse>=rab))/rnde; 

        itr_flagrab=rnde;

        % sub-step 2, max.perm
        if rnd>rnde 
            if pvrab<=thr
                disp('incremental permutation begins......');

                i=rnde+1;
                rabs=nan(1,rnd); 
                rabs(1:i-1)=rabse; 
                clear rabse;
                
                if isequal(cons,all)==1 
                    if sum(sum(isnan(dset)))==0

                       nant=nan(1,rnd-rnde);
                       rabs(rnde+1:end)=nant;
                       clear nant;
                            
                       switch rab2 % common cluster validation
                           case 1
                               rabs(rnde+1:end)=prod(dset(:,rnde+1:end))./(one.*ones(1,rnd-rnde)); 
                           case 2
                               rabs(rnde+1:end)=sum(dset(:,rnde+1:end))./(one.*ones(1,rnd-rnde));
                           case 3 
                               rabs(rnde+1:end)=mean(dset(:,rnde+1:end))./(one.*ones(1,rnd-rnde));
                           case 4 
                               rabs(rnde+1:end)=geomean(dset(:,rnde+1:end))./(one.*ones(1,rnd-rnde));
                           case 5
                               rabs(rnde+1:end)=conslen*(conslen-1)*prod(dset(:,rnde+1:end))./(one.*ones(1,rnd-rnde)); 
                           case 6
                               det=dset(:,rnde+1:end); 
                               det(det==0)=sv;
                               rabs1=zeros(1,rnd-rnde); 
                               if isequal(conslen,cond)==0
                                   error('data sets G+ size mismatch');
                               end                             
                               for k=1:conslen-1 
                                   for j=k+1:conslen 
                                       rabs1=dset(cons(k),rnde+1:end)./det(cons(j),:)+dset(cons(j),rnde+1:end)./det(cons(k),:)+rabs1;
                                   end
                               end                                 
                               rabs(rnde+1:end)=conslen*(conslen-1)*rabs1./(one.*ones(1,rnd-rnde));
                               clear rabs1 det;
                           case 7
                               det=dset(:,rnde+1:end);
                               det(det==0)=sv;
                               rabs1=zeros(1,rnd-rnde);
                               if isequal(conslen,cond)==0
                                   error('data sets G+ size mismatch');
                               end                             
                               for k=1:conslen-1
                                   for j=k+1:conslen 
                                       rabs1=dset(cons(k),rnde+1:end)./det(cons(j),:)+dset(cons(j),rnde+1:end)./det(cons(k),:)+rabs1;
                                   end
                               end
                               rabs1(rabs1==0)=sv;
                               rabs(rnde+1:end)=(sum(dset(cons,rnde+1:end))./rabs1)./(one./one);
                               clear rabs1 det;
                           otherwise
                               error('Invalid input (1-7) for cluster quality');
                       end

                       if sum(isnan(rabs))>0 || sum(isinf(rabs))>0
                            error('mistake in incremental permutation test');
                       end

                       pvrab=length(find(rabs>=rab))/rnd;
                       itr_flagrab=rnd;

                    else
                        error('Mistake: permutation test for the common clsuter');
                    end
        
                else 
                    if conslen>=cond || stoplen<=0
                        error('data size mismatch');
                    end
                    if cond~=(conslen+stoplen) || conslen<=0
                        error('data size mismatch');
                    end
                    
                    % full permutation in sub-step 1
                    if sum(sum(isnan(dset)))==0
                        % for differential cluster
                        switch rab2
                            case 1
                                det=dset(constop,rnde+1:end);
                                det(det==0)=sv;
                                if stoplen>1 
                                    det=prod(det);
                                end 
                                if conslen>1
                                    rabs(rnde+1:end)=prod(dset(cons,rnde+1:end))./det;
                                else
                                    rabs(rnde+1:end)=dset(cons,rnde+1:end)./det;
                                end
                                clear det;
                            case 2
                                if stoplen==1 
                                    det=dset(constop,rnde+1:end); 
                                 else 
                                    det=sum(dset(constop,rnde+1:end));
                                end
                                det(det==0)=sv;
                                if conslen>1
                                    rabs(rnde+1:end)=sum(dset(cons,rnde+1:end))./det;
                                else 
                                    rabs(rnde+1:end)=dset(cons,rnde+1:end)./det;
                                end 
                                clear det;
                            case 3
                                if stoplen==1 
                                    det=dset(constop,rnde+1:end); 
                                 else 
                                    det=mean(dset(constop,rnde+1:end));
                                end
                                det(det==0)=sv;
                                if conslen>1
                                    rabs(rnde+1:end)=mean(dset(cons,rnde+1:end))./det;
                                else 
                                    rabs(rnde+1:end)=dset(cons,rnde+1:end)./det;
                                end 
                                clear det;                    
                            case 4
                                if stoplen==1 
                                    det=dset(constop,rnde+1:end);
                                 else 
                                    det=geomean(dset(constop,rnde+1:end));
                                end
                                det(det==0)=sv;
                                if conslen>1
                                    rabs(rnde+1:end)=geomean(dset(cons,rnde+1:end))./det;
                                else 
                                    rabs(rnde+1:end)=dset(cons,rnde+1:end)./det;
                                end 
                                clear det;
                            case 5
                                rabse1=nan(1,rnde);
                                rabse2=nan(1,rnde);
                                if conslen>1
                                    rabse1=conslen*(conslen-1)*prod(dset(cons,rnde+1:end));
                                else
                                    rabse1=dset(cons,rnde+1:end);
                                end
                                det=dset(constop,rnde+1:end);
                                det(det==0)=sv;
                                if stoplen>1
                                    rabse2=stoplen*(stoplen-1)*prod(det);
                                else 
                                    rabse2=det; 
                                end
                                rabs(rnde+1:end)=rabse1./rabse2;
                                clear rabse1 rabse2 det;        
                            case 6
                                det=dset(:,rnde+1:end); 
                                det(det==0)=sv;
                                rabse1=nan(1,rnd-rnde);
                                if conslen>1
                                    rabse1=zeros(1,rnd-rnde);
                                    for k=1:conslen-1
                                       for j=k+1:conslen 
                                           rabse1=dset(cons(k),rnde+1:end)./det(cons(j),:)+dset(cons(j),rnde+1:end)./det(cons(k),:)+rabse1;
                                       end
                                    end
                                    rabse1=conslen*(conslen-1)*rabse1;
                                else 
                                    rabse1=dset(cons,rnde+1:end);
                                end                            
                                rabse2=nan(1,rnd-rnde);
                                if stoplen>1 
                                    rabse2=zeros(1,rnd-rnde);
                                    for k=1:stoplen-1
                                       for j=k+1:stoplen
                                           rabse2=dset(constop(k),rnde+1:end)./det(constop(j),:)+dset(constop(j),rnde+1:end)./det(constop(k),:)+rabse2;
                                       end
                                    end
                                    rabse2(rabse2==0)=sv; 
                                    rabse2=stoplen*(stoplen-1)*rabse2;
                                else 
                                    rabse2=det(constop,:); 
                                end
                                rabs(rnde+1:end)=rabse1./rabse2;
                                clear rabse1 rabse2 det;
                            case 7
                                det=dset(:,rnde+1:end);
                                det(det==0)=sv;
                                rabse1=nan(1,rnd-rnde);
                                if conslen>1
                                    det1=sum(dset(cons,rnde+1:end));
                                    rabse1=zeros(1,rnd-rnde);
                                    for k=1:conslen-1
                                       for j=k+1:conslen 
                                           rabse1=dset(cons(k),rnde+1:end)./det(cons(j),:)+dset(cons(j),rnde+1:end)./det(cons(k),:)+rabse1;
                                       end
                                    end
                                    rabse1(rabse1==0)=sv;
                                else 
                                    det1=dset(cons,rnde+1:end);
                                    rabse1=det(cons,:); 
                                end
                                rabse1=det1./rabse1;
                                rabse2=nan(1,rnd-rnde);
                                if stoplen>1  
                                    det2=sum(dset(constop,rnde+1:end));
                                    det2(det2==0)=sv;
                                    rabse2=zeros(1,rnd-rnde);
                                    for k=1:stoplen-1
                                       for j=k+1:stoplen
                                           rabse2=dset(constop(k),rnde+1:end)./det(constop(j),:)+dset(constop(j),rnde+1:end)./det(constop(k),:)+rabse2;
                                       end
                                    end
                                else
                                    det2=det(constop,:); 
                                    rabse2=dset(constop,rnde+1:end);
                                end
                                rabse2=det2./rabse2;
                                rabs(rnde+1:end)=rabse1./rabse2;
                                clear rabse1 rabse2 det det1 det2;
                            otherwise
                                error('Invalid input (1-7) for cluster quality');
                        end
                        
                        if sum(isnan(rabs))>0 || sum(isinf(rabs))>0
                            error('mistake in incremental permutation test');
                        end
                        pvrab=length(find(rabs>=rab))/rnd;
                        itr_flagrab=rnd;
                    
                    else
                        
                        if i~=rnde+1
                            error('loop index i: value inconsistent during the loop');
                        end

                        dst=nan(1,cond);
                        while i<=rnd               
                            pr=randperm(N);
                            for j=1:cond
                                Ds{j}=Dc{j}(pr,pr); 
                            end

                            if ctrig==2 
                                di=nan(1,cond);
                                do=nan(1,cond);
                                for j=1:cond
                                    di(j)=sum(sum(Ds{j}(1:intv,1:intv),2));
                                    do(j)=(sum(sum(Ds{j},2))-di(j))./(2*(N-intv));
                                end
                                do(do==0)=sv;
                                di=di./(2*intv);
                                dst=di./do;
                                if sum(isnan(dst))>0
                                     error('nan in dst');
                                end       
                            else
                                for j=1:cond
                                    dst(j)=mean(mean(Ds{j}(1:intv,1:intv),2));
                                end
                            end
                            switch rab2 
                                case 1
                                    rabs(i)=prod(dst(cons))./prod(dst(constop)); 
                                    if isinf(rabs(i))==1 
                                        det=dst(constop);
                                        det(det==0)=sv; 
                                        rabs(i)=prod(dst(cons))./prod(det);
                                    end
                                    clear det;
                                case 2
                                    if sum(dst(constop))==0
                                        rabs(i)=sum(dst(cons))./sv;
                                    else
                                        rabs(i)=sum(dst(cons))./sum(dst(constop));
                                    end
                                case 3
                                    if mean(dst(constop))==0
                                        rabs(i)=mean(dst(cons))./sv;
                                    else
                                        rabs(i)=mean(dst(cons))./mean(dst(constop));
                                    end 
                                case 4
                                    if geomean(dst(constop))==0
                                        rabs(i)=geomean(dst(cons))./sv;
                                    else
                                        rabs(i)=geomean(dst(cons))./geomean(dst(constop));
                                    end 
                                case 5
                                    if conslen>1
                                        rab1s=conslen*(conslen-1)*prod(dst(cons));
                                    else 
                                        rab1s=prod(dst(cons)); 
                                    end
                                    det=dst(constop);
                                    det(det==0)=sv;
                                    if stoplen>1
                                        rab2s=stoplen*(stoplen-1)*prod(det);
                                    else 
                                        rab2s=prod(det); 
                                    end
                                    rabs(i)=rab1s./rab2s;
                                    clear rab1s rab2s det;
                                case 6
                                    det=dst;
                                    det(det==0)=sv;
                                    if conslen>1  
                                        rab1s=0;
                                        for j=1:conslen-1
                                           for k=j+1:conslen
                                               rab1s=dst(cons(j))./det(cons(k))+dst(cons(k))./det(cons(j))+rab1s;
                                           end
                                        end
                                        rab1s=conslen*(conslen-1)*rab1s;
                                    else 
                                        rab1s=dst(cons);
                                    end                      
                                    if stoplen>1
                                        rab2s=0;
                                        for j=1:stoplen-1
                                           for k=j+1:stoplen
                                               rab2s=dst(constop(j))./det(constop(k))+dst(constop(k))./det(constop(j))+rab2s;
                                           end
                                        end
                                        if rab2s==0
                                            rab2s=sv;
                                        end
                                        rab2s=stoplen*(stoplen-1)*rab2s;
                                    else 
                                        rab2s=det(constop); 
                                    end
                                    rabs(i)=rab1s./rab2s;
                                    clear rab1s rab2s det;
                                case 7
                                    dets=dst;
                                    dets(dets==0)=sv;
                                    if conslen>1 
                                        rab1s=0;
                                        for j=1:conslen-1
                                           for k=j+1:conslen
                                               rab1s=dst(cons(j))./dets(cons(k))+dst(cons(k))./dets(cons(j))+rab1s;
                                           end
                                        end
                                        if rab1s==0 
                                            rab1s=sv;
                                        end
                                    else 
                                        rab1s=dets(cons); 
                                    end
                                    rab1s=sum(dst(cons))./rab1s; 
                                    if stoplen>1
                                        rab2s=0;
                                        for j=1:stoplen-1
                                           for k=j+1:stoplen
                                               rab2s=dst(constop(j))./dets(constop(k))+dst(constop(k))./dets(constop(j))+rab2s;
                                           end
                                        end
                                    else 
                                        rab2s=dst(constop); 
                                    end
                                    det=sum(ds(constop)); 
                                    det(det==0)=sv;
                                    rab2s=det./rab2s; 
                                    rabs(i)=rab1s./rab2s; 
                                    clear rab1s rab2s det dets;
                                otherwise
                                    error('Invalid input (1-7) for cluster quality');
                            end 

                            if isnan(rabs(i))==1 || isinf(rabs(i))==1
                                error('mistake in Validation sub-step 2 incremental permutation test');
                            end

                            pvrab=length(find(rabs(1:i)>=rab))/i; 
                            itr_flagrab=i;
                            if pvrab>thr
                                i=rnd+1; 
                            else 
                                i=i+1; 
                            end
                        end 
                        
                        if itr_flagrab<rnd
                            if length(rabs)==rnd
                                rabs(itr_flagrab+1:end)=[]; 
                            else
                                error('mistake: initializing rabs');
                            end
                        end 
                    end
                end

            else 
                disp('the cluster found in Validation sub-step 1 is not significant (based on Validation sub-step 2 (min.perm))');
            end
        else
            error('min.perm can not be higher than max.perm');
        end

    else     
        disp('the cluster presents in no conditions');
        rab=nan(1); 
        pvrab=nan(1); 
        itr_flagrab=nan(1); 
        consout=cons; 
    end
    
    disp('Validation Step 2 done ....');
    fprintf('the overall p-value is %f for cluster (node %d to %d, %d permutations\n',pvrab,td1,td2,itr_flagrab);
    
end

end % end of validation









    
    
    
    
    
    
    













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

end % end of C3Dfunction [V,Lamda]=e_trans(D)
% Compute S and the corresponding eigensystems for the symmetric(co-expression) input data matrices D
% Input:
%     D:    each cell D{i} is a symmetric co-expression matrix
% Output:
%     V:    the common factor of hosgvd decomposition of D, 
%     each column is an eigenvector   
%     Lamda: each diagonal value is an eignevalue 
% W is computed directly from the D
%
% Written by Xiaolin Xiao, 2012


if nargin < 1 || isempty(D)
    error(message('TooFewInputs'));
end

cond=length(D);
if cond<2
    error('number of conditions is less than 2');
end

drsize=nan(1,cond);
real_tag=nan(1,cond);
for i=1:cond
    drsize(i)=size(D{i},1);
    real_tag(i)=isreal(D{i});
end
if isequal(drsize, drsize(1).*ones(1,cond))==0
    error(message('MATLAB:hogsvd:MatrixRowMismatch'))
end
if isequal(real_tag,ones(1,cond))==0
    error('Complex numbers exist');
end


A=D; % take co-expression

S=zeros(length(A{1}));
for i=1:cond-1
    for j=i+1:cond
        if size(D{j},1)<size(D{j},2)
            if size(D{i},1)<size(D{i},2)
                S=S+(A{i}*fast_pinv(D{j}')*fast_pinv(D{j})+A{j}*fast_pinv(D{i}')*fast_pinv(D{i}));
            else
                S=S+(A{i}*fast_pinv(D{j}')*fast_pinv(D{j})+A{j}*pinv(D{i}')*pinv(D{i}));
            end
        else
            if size(D{i},1)<size(D{i},2)
                S=S+(A{i}*pinv(D{j}')*pinv(D{j})+A{j}*fast_pinv(D{i}')*fast_pinv(D{i}));
            else
                S=S+(A{i}*pinv(D{j}')*pinv(D{j})+A{j}*pinv(D{i}')*pinv(D{i}));
            end
        end           
    end
end
S=S./(cond*(cond-1));
rs=rank(S);

if rs<length(S)-1
    [V,Lamda]=eigs(S,rs);
else
    [V,Lamda]=eig(S);
end

if isreal(Lamda)==0
    lam=diag(Lamda);
    k=1;
    while isreal(lam(k))==1
        k=k+1;
    end
    if k<=1
        error('non real patterns found');
    else
        Lamda=Lamda(1:k-1,1:k-1);
        V=V(:,1:k-1);
    end
end


colnorm=sqrt(sum(V.^2,1)); 
colnorm=repmat(colnorm,[size(V,1) 1]);
V=V./colnorm; 


end % end of function








function [V,Lamda]=hogsvd_trans(D)
% Compute the eigensystems under the Higher Order GSVD for general real multiple data matrices
% Input:
%     D:    each cell D{i} is one condition of the real multiple input data
%     matrices (each cell must be a real-valued data matrix with the same 
%     ROW size)
% Output:
%     V:    the common factor of hogsvd decomposition of D, each
%     column is an eigenvector
%     Lamda:each diagonal value is an eignevalue 
%
% Written by Xiaolin Xiao, 2012


if nargin < 1 || isempty(D)
    error(message('TooFewInputs'));
end

% cell numbers: number of matrices
cond=length(D);
if cond<2
    error('number of conditions is less than 2');
end

drsize=nan(1,cond);
real_tag=nan(1,cond);
for i=1:cond
    drsize(i)=size(D{i},1);
    real_tag(i)=isreal(D{i});
end
if isequal(drsize, drsize(1).*ones(1,cond))==0
    error(message('MATLAB:hogsvd:MatrixRowMismatch'))
end
if isequal(real_tag,ones(1,cond))==0
    error('Complex numbers exist');
end

A=cell(cond,1);
for i=1:cond
    A{i}=D{i}*D{i}'; 
end

S=zeros(length(A{1}));
for i=1:cond-1
    for j=i+1:cond
        if size(D{j},1)<size(D{j},2)
            if size(D{i},1)<size(D{i},2)
                S=S+(A{i}*fast_pinv(D{j}')*fast_pinv(D{j})+A{j}*fast_pinv(D{i}')*fast_pinv(D{i}));
            else
                S=S+(A{i}*fast_pinv(D{j}')*fast_pinv(D{j})+A{j}*pinv(D{i}')*pinv(D{i}));
            end
        else
            if size(D{i},1)<size(D{i},2)
                S=S+(A{i}*pinv(D{j}')*pinv(D{j})+A{j}*fast_pinv(D{i}')*fast_pinv(D{i}));
            else
                S=S+(A{i}*pinv(D{j}')*pinv(D{j})+A{j}*pinv(D{i}')*pinv(D{i}));
            end
        end           
    end
end
S=S./(cond*(cond-1));
rs=rank(S);

if rs<length(S)-1
    [V,Lamda]=eigs(S,rs);
else
    [V,Lamda]=eig(S);
end

if isreal(Lamda)==0
    lam=diag(Lamda);
    k=1;
    while isreal(lam(k))==1
        k=k+1;
    end
    if k<=1
        error('non real patterns found');
    else
        Lamda=Lamda(1:k-1,1:k-1);
        V=V(:,1:k-1);
    end
end


colnorm=sqrt(sum(V.^2,1)); 
colnorm=repmat(colnorm,[size(V,1) 1]);
V=V./colnorm; 


end % end of function hogsvd








function Y = fast_pinv(G)
% Fast computation of moore-penrose inverse of matrix G(m*n) if m<n
% reference: P. Courrieu, Fast Computation of Moore-Penrose Inverse
% Matrices, Neural Information Processing -Letters and Reviews, 8 (2),
% 2005.


[m,n]=size(G); 
transpose=false;
if m<n
    transpose=true;
    A=G*G';
    n=m;
else
    A=G'*G;
end

dA=diag(A); 
tol= min(dA(dA>0))*1e-9;
L=zeros(size(A));
r=0;
for k=1:n
    r=r+1;
    L(k:n,r)=A(k:n,k)-L(k:n,1:(r-1))*L(k,1:(r-1))';
    if L(k,r)>tol
        L(k,r)=sqrt(L(k,r));
        if k<n
            L((k+1):n,r)=L((k+1):n,r)/L(k,r);
        end
    else
        r=r-1;
    end
end
L=L(:,1:r);

M=inv(L'*L);
if transpose
    Y=G'*L*M*M*L';
else
    Y=L*M*M*L'*G';
end#! /usr/bin/Rscript
# written by Xiaolin Xiao, 2013 

getwd()
a=getwd()
library(fdrtool)
fin=paste(a,"vx.txt",sep="/")
vin=read.table(fin)
vin=vin[,1]
fdr.out=fdrtool(vin,statistic="normal",plot=FALSE)
qval_all=fdr.out$qval
fout_qall=paste(a,"vx_outqval_all.txt",sep="/")
write.table(qval_all,fout_qall)

