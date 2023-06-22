

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

