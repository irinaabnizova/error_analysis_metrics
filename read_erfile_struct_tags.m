function [maQ,ncyc,coH1,coH2,coL1,coL2,winHR1,winHR2,winLR1,winLR2,msh_1,msh_2,msl_1,msl_2,E12,set_1,set_2,rch,rcl,mimh,miml,prcnh,prcnl]=read_erfile_struct_tags(ef,tags)
%============== if error file exists, read it into arrays 
%and structures : rch,
  
%-----------------------initiate values of arrays
maQ=0;
ncyc=0;
coH1=[];
coH2=[];
coL1=[];
coL2=[];
winLR1=[];
winLR2=[];
winHR1=[];
winHR2=[];
E12=[];

msh_1=[];
msh_2=[];
msl_1=[];
msl_2=[];
set_1=[];
set_2=[];

rch=[];rcl=[];mimh=[];miml=[];prcnh=[];prcnl=[];

% for all cases in tags=names from error file
i=0;
%ef
if ef > 0,

  %----------------------------tag1=SET

  i=i+1;
  if i==1,
  tag=tags(1);
  [sub,read_Q,maQ]=read_set(ef,tag);
  siSET=size(sub);
  si_RQ=size(read_Q);

  ma=max(siSET);
    if ma>0,
    set_1=sub(1:ma/2,:);
    set_2=sub(ma/2+1:ma,:);
    totEa_R1=sum(sum(set_1));
    totEa_R2=sum(sum(set_2));
    E12=[totEa_R1,totEa_R2];
    end
  end% i==1


%2----------------read the parts of error file about RCH RCL (gen sub counts HQ)
 i=i+1;
 frewind(ef);
 if i==2,
    tag=tags(i);
    tag_current=tags(i);
    [msH,ncyc]=tag_read_rc(ef,tag,tag_current);
    si_subH=size(msH);
    if max(si_subH)>0,
    coH1=msH(ncyc+1,:);
    coH2=msH(2*(ncyc+1),:);
    msh_1=msH(1:ncyc,:);
    msh_2=msH(ncyc+2:2*ncyc+1,:);
    end
    %=======================rch struct
    field1='read'; value1={1,2};
    field2='count'; value2={coH1 coH2};
    rch=struct(field1,value1,field2,value2);
    
    %=======mimh struct
   field1='read'; value1={1,2};
   field3='count'; value3={msh_1 msh_2};
   field2='cycle'; value2={1:ncyc};
   mimh=struct(field1,value1,field2,value2,field3,value3);


 end

 %-----------------tag3=RCL
 i=i+1;
 frewind(ef);
 if i==3,
    tag=tags(i);
    tag_current=tags(i);
    [msL,ncyc]=tag_read_rc(ef,tag,tag_current);
    si_subL=size(msL);
    if max(si_subL)>0,
    coL1=msL(ncyc+1,:);
    coL2=msL(2*(ncyc+1),:);
    msl_1=msL(1:ncyc,:);
    msl_2=msL(ncyc+2:2*ncyc+1,:);
    end
    
    %=======================rch struct
    field1='read'; value1={1,2};
    field2='count'; value2={coL1 coL2};
    rcl=struct(field1,value1,field2,value2);
  
   %=======miml struct
   field1='read'; value1={1,2};
   field3='count'; value3={msl_1 msl_2};
   field2='cycle'; value2={1:ncyc};
   miml=struct(field1,value1,field2,value2,field3,value3);
   
 end

 %--------------------------------tag4= PRCNH
 i=i+1;
 frewind(ef);
 if i==4,
    tag=tags(i);
    tag_current=tags(i);
    [subwH,len_r]=tag_read_win(ef,tag,tag_current);
    si_subwH=size(subwH);

    if max(si_subwH)>0,
    winHR11=subwH(1:len_r,:);
    winHR21=subwH(len_r+1:2*len_r,:);
    [winHR1]=stripped_convert_new(winHR11);
    [winHR2]=stripped_convert_new(winHR21);
    end
    
        
    %=======================prcnh struct
    field1='read'; value1={1,2};
    field2='count'; value2={winHR1 winHR2};
    prcnh=struct(field1,value1,field2,value2);
    
 end
 %--------------------------------tag5= PRCNL
 i=i+1;
 frewind(ef);
 if i==5,
    tag=tags(i);
    tag_current=tags(i);
    [subwL,len_r]=tag_read_win(ef,tag,tag_current);
    si_subwL=size(subwL);

    if max(si_subwL)>0,
    winLR11=subwL(1:len_r,:);
    winLR21=subwL(len_r+1:2*len_r,:);
    [winLR1]=stripped_convert_new(winLR11);
    [winLR2]=stripped_convert_new(winLR21);
    end
    
    %=======================prcnh struct
    field1='read'; value1={1,2};
    field2='count'; value2={winLR1 winLR2};
    prcnl=struct(field1,value1,field2,value2);
    
 end
 
 
else
    display(' error message: no error file, exit with initiated empty arrays values');
    return
end
end% function

%==================================define subfunctions
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++1
function [sub,read_Q,maQ]=read_set(ef,tag)
% read 'SET' part from error file

sub=[]; 
read_Q=[];
maQ=0;

r=[];
q=[];

tf = strcmp('SET',tag);
if (tf ) < 1,
    display('wrong tag, not SET ');
    return
else 
    
i=1;
k=0;
while ~feof(ef),
line_ex = fgetl(ef);
  if line_ex(1)=='#',
      a=1;
  else
      C = strsplit(line_ex);
      tf = strcmp(C(1),tag);
      if tf > 0,
          k=k+1;
         b = str2double(C(2));
         c = str2double(C(3));    
         r(k)=b;
         q(k)=c;
         xj=5:2:length(C);
      
         for i1=1:length(xj),
             su(i1)=str2double(C(xj(i1)));
         end
            sub(k,:)=su;
      end
      
  end
  i=i+1;
end

m=0;
for i=1:length(r),
    if r(i)==1.0,
        m=m+1;
        qq(m)=q(i);
    end
end
maQ=max(q);

%----------------fill in new array of qualities
read_Q =[r' q'];
end       
end
 
%=================================================2
function [sub,ncyc]=tag_read_rc(ef,tag,tag_current);
% reads RCH and RCL parts of error file; outputs array (12x16) and structure RC

sub=[];
ncyc=0;

tf = strcmp(tag_current,tag);
if (tf ) <1,
    display('wrong tag, not current tag ');
    return;
else 
    
i=1;
k=0;
while ~feof(ef),
line_ex = fgetl(ef);
  if line_ex(1)=='#',
  else
      C = strsplit(line_ex);
      tf = strcmp(C(1),tag);
      if tf > 0,
          k=k+1;
         b = str2double(C(2));
         c = str2double(C(3));    
         r(k)=b;
         cyc(k)=c;
         xj=5:2:length(C);
       
         for i1=1:length(xj),
             su(i1)=str2double(C(xj(i1)));
         end
             sub(k,:)=su;
      end
      
  end
  i=i+1;
end

if k>0,
 m=0;
 for i=1:length(r),
    if r(i)==1.0,
        m=m+1;
        ccyc(m)=cyc(i);
    end
 end
 ncyc=max(ccyc)+1;

end
end 
end
%============================================3
function [sub,len_r]=tag_read_win(ef,tag,tag_current);
% reads the part of error file with PRCNH and PRCNL;outputs 16x12 array 

sub=[];
len_r=0;
r=[];

tf = strcmp(tag_current,tag);
if (tf ) <1,
    display('wrong tag, not current tag ');
    return;
else 
    
    
i=1;
k=0;
while ~feof(ef),
line_ex = fgetl(ef);
  if line_ex(1)=='#',
  else
      C = strsplit(line_ex);
      tf = strcmp(C(1),tag);
      if tf > 0,
          k=k+1;
         b = str2double(C(2));
         c = str2double(C(3));    
         r(k)=b;
         cyc(k)=c;
         xj=4:2:length(C);
       
         for i1=1:length(xj),
             su(i1)=str2double(C(xj(i1)));
         end
         
         sub(k,:)=su;
      end
      
  end
  i=i+1;
end

m=0;
if size(r)>0,
  for i=1:length(r),
    if r(i)==1.0,
        m=m+1;
        
    end
  end
end
len_r=m;

end  
end
   %=======================================5
function [new_win]=stripped_convert_new(win_R1);
%coverts win_R1 (as in error quality) into converient format, new_win,
%where each of 12 columns are counts of windows X.Y for each of 12 substitutions

%INPUT input win_R1 (16x12) from error file  see below:

%# Effect of previous base and next base high predictor. Use `grep ^PNCRH | cut -f 2-` to extract this part
%# 16 rows: one row per pair= previous/next base, columns (12 for each subs) is previous base+substitution+next base and count for 12 substitutions

%PCRNH	1  caA bw AX	AACA	6	AACC	10	AACG	6	AACT	2	AAGA	159	AAGC	125	AAGG	133	AAGT	102	AATA	45	AATC	51	AATG	52	AATT	65
%PCRNH	1  caC	        ACAA	11	ACAC	23	ACAG	10	ACAT	12	ACGA	2	ACGC	12	ACGG	8	ACGT	4	ACTA	129	ACTC	156	ACTG	276	ACTT	74
%PCRNH	1  caG	        AGAA	51	AGAC	179	AGAG	77	AGAT	53	AGCA	45	AGCC	32	AGCG	20	AGCT	13	AGTA	646	AGTC	701	AGTG	398	AGTT	209
%PCRNH	1  caT	        ATAA	49	ATAC	64	ATAG	49	ATAT	53	ATCA	116	ATCC	92	ATCG	137	ATCT	119	ATGA	2	ATGC	7	ATGG	8	ATGT	5

%PCRNH	1  caA bw CX	CACA	6	CACC	4	CACG	14	CACT	12	CAGA	146	CAGC	65	CAGG	155	CAGT	124	CATA	60	CATC	39	CATG	60	CATT	72
%PCRNH	1	CCAA	56	CCAC	25	CCAG	22	CCAT	30	CCGA	13	CCGC	3	CCGG	7	CCGT	7	CCTA	175	CCTC	34	CCTG	138	CCTT	66
%PCRNH	1	CGAA	129	CGAC	363	CGAG	97	CGAT	244	CGCA	33	CGCC	37	CGCG	36	CGCT	22	CGTA	463	CGTC	509	CGTG	306	CGTT	251
%PCRNH	1	CTAA	23	CTAC	46	CTAG	37	CTAT	45	CTCA	178	CTCC	118	CTCG	204	CTCT	202	CTGA	3	CTGC	2	CTGG	9	CTGT	10

%PCRNH	1  caA bw GX	GACA	3	GACC	7	GACG	6	GACT	4	GAGA	99	GAGC	78	GAGG	94	GAGT	83	GATA	60	GATC	58	GATG	46	GATT	72
%PCRNH	1	GCAA	97	GCAC	20	GCAG	20	GCAT	29	GCGA	11	GCGC	11	GCGG	10	GCGT	14	GCTA	201	GCTC	192	GCTG	307	GCTT	162
%PCRNH	1	GGAA	53	GGAC	214	GGAG	36	GGAT	94	GGCA	33	GGCC	21	GGCG	7	GGCT	13	GGTA	348	GGTC	388	GGTG	145	GGTT	184
%PCRNH	1	GTAA	26	GTAC	47	GTAG	33	GTAT	51	GTCA	3520	GTCC	96	GTCG	71	GTCT	104	GTGA	4	GTGC	5	GTGG	11	GTGT	22

%PCRNH	1 caA bw TX	TACA	4	TACC	3	TACG	2	TACT	3	TAGA	96	TAGC	4020	TAGG	209	TAGT	90	TATA	63	TATC	46	TATG	16	TATT	49
%PCRNH	1	TCAA	28	TCAC	11	TCAG	19	TCAT	23	TCGA	15	TCGC	5	TCGG	7	TCGT	5	TCTA	186	TCTC	65	TCTG	141	TCTT	84
%PCRNH	1	TGAA	134	TGAC	272	TGAG	112	TGAT	86	TGCA	74	TGCC	47	TGCG	46	TGCT	14	TGTA	774	TGTC	822	TGTG	512	TGTT	211
%PCRNH	1	TTAA	45	TTAC	63	TTAG	31	TTAT	43

new_win=[];
N=4;

a=1;
for f=1:3,
 for k=1:4,
      j=a+(k-1)*N;
      for t1=1:4,
            new_win(j+t1-a,f)=win_R1(j,t1+(f-1)*N); 
      end
  end
end 


c=2;
for f=1:3,
  for k=1:4,
      j=c+(k-1)*N;
      for t1=1:4,
            new_win(j+t1-c,f+2*c-1)=win_R1(j,t1+(f-1)*N); 
      end
  end
end 

g=3;
for f=1:3,
   for k=1:4,
      j=g+(k-1)*N;
      for t=1:4,
            new_win(j+t-g,f+2*g)=win_R1(j,t+(f-1)*N);
      end
   end
end 


t=4;
for f=1:3,
     for k=1:4,
      j=t+(k-1)*N;% caG  3,
      for t1=1:4,
            new_win(j+t1-t,f+2*t+1)=win_R1(j,t1+(f-1)*N); 
      end
     end
 end 
end% subfunction


        
        
        
        
        
        
        
        
        
