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









    
    
    
    
    
    
    













