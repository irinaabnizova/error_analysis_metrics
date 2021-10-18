   function [name_out]=compute_visualise_percycle(file_in,name,thr_infla)
   
       thr_infla

       tags={'SET','RCH','RCL','PRCNH','PRCNL'};

       %name='38912/38912_4#13'
       s2 = regexp(name, '/', 'split');
       name1=s2{2}
       %qn = sprintf('%s_quality_error.txt',name);= file_in
       %qn = sprintf('%s_substitution_analysis.txt',name);
       ef=fopen(file_in);
    
      if ef>0,
    %2.1 ============== if error file exists, read it into arrays 
                    %and structures : rch,rcl,mimh,miml,prcnh,prcnl

          [maQ,ncyc,coH1,coH2,coL1,coL2,winHR1,winHR2,winLR1,winLR2,msh_1,msh_2,msl_1,msl_2,E12]=read_erfile_struct_tags(ef,tags);
     
    %2.2 ---------------------------GP plots
          [coH1p,coL1p,coH2p,coL2p,mahh,mall,casm_plot,total]=general_profile_plots(coL1,coH1,coL2,coH2,name1);
        
    %3. ================compute metrics with correction of Ti asymmymetry
    %if needed- SKIP here
                  
                  
    %3. PER CYCLE -------------------------per cycle metrics and plots
         [mshp_1,mshp_2, mslp_1, mslp_2,maL1,maH1]=metrics_per_cycle_plots_R1R2_HL_trim(msh_1,msh_2,msl_1,msl_2);
         title([' Per cycle mismatch,  HQ R1, ',name1],'Interpreter','none','FontSize',14);
 
           %-------------------check for inflated Ti HQ subs per cycle zones
                      
           %3.2.1-------------first 10,20,30 cycles
            subs_num=subs_numbers;        
            ind_Ti=[2,6,7,11];
            sub_Ti=subs_num(ind_Ti);
           
           j=ind_Ti(1); %AG;   %j=6 % CT 7-GA 11 TC
           [ratio_fs_1,ratio_ft_1,ratio_fs_2,ratio_ft_2,inflaAG]=inflated_first_cyclesLO(thr_infla,mshp_1,mshp_2,j,subs_num);
           j=ind_Ti(2);;%CT
           [ratio_fs_1,ratio_ft_1,ratio_fs_2,ratio_ft_2,inflaCT]=inflated_first_cyclesLO(thr_infla,mshp_1,mshp_2,j,subs_num);
           j=ind_Ti(3);;%GA
           [ratio_fs_1,ratio_ft_1,ratio_fs_2,ratio_ft_2,inflaGA]=inflated_first_cyclesLO(thr_infla,mshp_1,mshp_2,j,subs_num);
           j=ind_Ti(4);;%TC
           [ratio_fs_1,ratio_ft_1,ratio_fs_2,ratio_ft_2,inflaTC]=inflated_first_cyclesLO(thr_infla,mshp_1,mshp_2,j,subs_num);
         
           infla_Ti_R12=[sub_Ti' [inflaAG;inflaCT;inflaGA;inflaTC]]
           
           %--------------------save
           [name_out]=save_cycle_metrics_tag(name,infla_Ti_R12);
    
    %=================================SUBFUNCTIONS
    else
        display('no such a data file');
    end % if empty
   
 end
   
   
   %==============================subfunctions
function [subs]=subs_numbers

%unfold in vector
subs(1)=12;% sub_matr(2,1)=slopp(1);
subs(2)=13;%sub_matr(3,1)=slopp(2);
subs(3)=14;%sub_matr(4,1)=slopp(3);%

subs(4)=21;%sub_matr(1,1)=slopp(4);
subs(5)=23;%sub_matr(3,2)=slopp(5);
subs(6)=24;%sub_matr(4,2)=slopp(6);

subs(7)=31;%sub_matr(1,2)=slopp(7);
subs(8)=32;%sub_matr(2,2)=slopp(8);
subs(9)=34;%sub_matr(4,3)=slopp(9);

subs(10)=41;%sub_matr(1,3)=slopp(10);
subs(11)=42;%sub_matr(2,3)=slopp(11);
subs(12)=43;%sub_matr(3,3)=slopp(12);
 
end

%=======subs2
function [ratio_fs_1,ratio_ft_1,ratio_fs_2,ratio_ft_2,infla,sub]=inflated_first_cyclesLO(thr_infla,mshp_1,mshp_2,j,subs)
 %3.2.1-------------first 10,20,30 cycles for one of 12 sub (j here)
 sub=subs(j);
 %-----------OUTPUT
 infla1=0;
 infla2=0;
 infla=[infla1,infla2];
 
 %ratio_fs_1,ratio_ft_1,ratio_fs_2,ratio_ft_2= ratios first/second 10
 %cycles
 
 %----------------input: subs j (2,6,7,11 are j for Ti)
           
           %------------------------AG 2
           first10_AG_1=sum(mshp_1(1:10,j));
           second10_AG_1=sum(mshp_1(11:20,j));
           third10_AG_1=sum(mshp_1(21:30,j));
           
           first10_AG_2=sum(mshp_2(1:10,j));
           second10_AG_2=sum(mshp_2(11:20,j));
           third10_AG_2=sum(mshp_2(21:30,j));
           
           fst_AG_1=[first10_AG_1,second10_AG_1,third10_AG_1];
           fst_AG_2=[first10_AG_2,second10_AG_2,third10_AG_2];
           
           ratio_fs_1=fst_AG_1(1)/fst_AG_1(2);
           ratio_ft_1=fst_AG_1(1)/fst_AG_1(3);
           
           %thr_infla=1.3; %thr inflation first cycles
           if  ratio_fs_1 > thr_infla & ratio_ft_1 > thr_infla,
               display('the first cycle Ti R1 are inflated');
               ratio_first_1=ratio_ft_1
               infla1=1;
           else
               %display('the first cycle Ti R1 are NOT inflated');
               infla1=0;
           end
           
           %--------------Read2
           ratio_fs_2=fst_AG_2(1)/fst_AG_2(2);
           ratio_ft_2=fst_AG_2(1)/fst_AG_2(3);
           
           %thr_infla=1.3 %thr inflation first cycles
           if  ratio_fs_2 > thr_infla & ratio_ft_2 > thr_infla,
               display('the first cycle Ti R2 are inflated');
               ratio_first_2=ratio_ft_2
               infla2=1;
           else
               %display('the first cycle Ti R2 are NOT inflated');
               infla2=0;
           end
           infla=[infla1,infla2];
end
  
%================subfunct
function [fn1]=save_cycle_metrics_tag(name,infla_Ti_R12)
%saves  Ti HQ cycle possible inflation metrics: 0 (No) or 1 (Yes), for
%Read1 (_R1) and Read2 (_R2)

%---------------------INPUT
% name = 38124/38124_2#154 (example, tag name with a path to directory name)

%infla_Ti_R12 =

%         13.00             0             0
%         24.00             0          1.00
%         31.00             0             0
%         42.00             0             0

%----------------OUTPUT
% name of an output file with a path: fn1 = sprintf('%s_Ti_first_cycle_metrics.txt',name);


   fn1 = sprintf('%s_Ti_first_cycle_metrics.txt',name);
   dense_cr3=fopen(fn1,'w');
   
   fprintf(dense_cr3,'#Ti HQ first cycle metrics for a given tag \n'); 
  
   fprintf(dense_cr3,'stats  AG_inflated_R1   %d  AG_inflated_R2   %d\n',infla_Ti_R12(1,2),infla_Ti_R12(1,3)); 
   fprintf(dense_cr3,'stats  CT_inflated_R1   %d  CT_inflated_R2   %d\n',infla_Ti_R12(2,2),infla_Ti_R12(2,3)); 
   fprintf(dense_cr3,'stats  GA_inflated_R1   %d  GA_inflated_R2   %d\n',infla_Ti_R12(3,2),infla_Ti_R12(3,3)); 
   fprintf(dense_cr3,'stats  TC_inflated_R1   %d  TC_inflated_R2   %d\n',infla_Ti_R12(4,2),infla_Ti_R12(4,3)); 
   
   fclose(dense_cr3);

 %F2===================make a file with main stats and metrics: numbers only

     fn = sprintf('%s_Ti_first_cycle_metrics_num.txt',name);
     dense_cr3=fopen(fn,'w');
   
   
    fprintf(dense_cr3,'#Ti HQ first cycle metrics for a given tag, num \n'); 
    fprintf(dense_cr3,'%d\t%d\t%d\n',infla_Ti_R12(1,1),infla_Ti_R12(1,2),infla_Ti_R12(1,3));
    fprintf(dense_cr3,'%d\t%d\t%d\n',infla_Ti_R12(2,1),infla_Ti_R12(2,2),infla_Ti_R12(2,3));
    fprintf(dense_cr3,'%d\t%d\t%d\n',infla_Ti_R12(3,1),infla_Ti_R12(3,2),infla_Ti_R12(3,3));
    fprintf(dense_cr3,'%d\t%d\t%d\n',infla_Ti_R12(4,1),infla_Ti_R12(4,2),infla_Ti_R12(4,3));
  
 
    fclose(dense_cr3);
end
%=======================subfuncts
function [coH1p,coL1p,coH2p,coL2p,mahh,mall,casm_plot,total]=general_profile_plots(coL1,coH1,coL2,coH2,name);
  coH1H2L1L2=[coH1' coH2' coL1' coL2'];
     total=sum(sum(coH1H2L1L2));
     maH=5;
     maL=16;
     %------------------Plot GEneral subs profile AND compute analog of CASM plot,plot it
     fact=100000/total;
    [coH1p,coL1p,coH2p,coL2p,mahh,mall,casm_plot]=plot_gen_prof_CASM(coL1,coH1,coL2,coH2,maH,maL,fact,name);
     title(['      Substitution profile, HQ R1, ',name],'Interpreter','none','FontSize',14);
 end
%=========================subfunctions
function [coH1,coL1,coH2,coL2,mah,mal,casm_plot]=plot_gen_prof_CASM(coL_1,coH_1,coL_2,coH_2,maH,maL,fact,name);

%-------------output 
%casm_plot=round([CA_chunk CG_chunk CT_chunk TA_chunk TC_chunk TG_chunk]*fact);

%------------------normalised by dividing on total bases:  normB

%-----------output  % of all freqencies
% compare two Reads general profiles:  coL_1,coH_1 and   coL_2,coH_2

%--------------------input : four 12-vectors of HQ and LQ substitution
%counts

%-----------------notations
z=1:12;
subs={'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG'};

subsCASM={'CA','CG','CT','TA','TC','TG'};
subsCASM_12={'CA1','CA2','GT1','GT2','CG1','CG2','GC1','GC2','CT1','CT2','GA1','GA2','TA1','TA2','AT1','AT2','TC1','TC2','AG1','AG2','TG1','TG2','AC1','AC2'};
zCAZM=1:24;
subsCASM_12_1={'CA1','CA2','CAr1','CAr2','CG1','CG2','CGr1','CGr2','CT1','CT2','CTr1','CTr2','TA1','TA2','TAr1','TAr2','TC1','TC2','TCr1','TCr2','TG1','TG2','TGr1','TGr2'};

%1------------------CASM plot: only HQ Q30+
ac1=coH_1(1);ac2=coH_2(1);
ag1=coH_1(2);ag2=coH_2(2);
at1=coH_1(3);at2=coH_2(3);

ca1=coH_1(4);ca2=coH_2(4);
cg1=coH_1(5);cg2=coH_2(5);
ct1=coH_1(6);ct2=coH_2(6);

ga1=coH_1(7);ga2=coH_2(7);
gc1=coH_1(8);gc2=coH_2(8);
gt1=coH_1(9);gt2=coH_2(9);

ta1=coH_1(10);ta2=coH_2(10);
tc1=coH_1(11);tc2=coH_2(11);
tg1=coH_1(12);tg2=coH_2(12);

%fact=0.001;
CA_chunk=[ca1,ca2,gt1,gt2];
CG_chunk=[cg1,cg2,gc1,gc2];
CT_chunk=[ct1,ct2,ga1,ga2];


TA_chunk=[ta1,ta2,at1,at2];
TC_chunk=[tc1,tc2,ag1,ag2];
TG_chunk=[tg1,tg2,ac1,ac2];



figure('units','normalized','outerposition',[0 0 1 1]);
bar([CA_chunk CG_chunk CT_chunk TA_chunk TC_chunk TG_chunk]*fact);
set(gca,'XTick',1:24);
set(gca,'XTickLabel',subsCASM_12_1,'FontSize',8, 'FontWeight','bold');
grid;
title(['      Substitution profile CASM-like, HQ, ',name],'Interpreter','none','FontSize',12);

casm_plot=round([CA_chunk CG_chunk CT_chunk TA_chunk TC_chunk TG_chunk]*fact);


%2=====================NPG-style general profile plot

tot1=sum(sum(coH_1))+sum(sum(coL_1));
tot2=sum(sum(coH_2))+sum(sum(coL_2));
coH1=100*coH_1/tot1;%(sum(coH_1)+sum(coL_1));
coL1=100*coL_1/tot1;%(sum(coH_1)+sum(coL_1));


coH2=100*coH_2/tot2;%(sum(coH_2)+sum(coL_2));
coL2=100*coL_2/tot2;%(sum(coH_2)+sum(coL_2));


mal=max(max(max([coL2;coL1])));
maL1=max(maL,mal);


mah=max(max(max([coH2;coH1])));
maH1=max(maH,mah);

%figure;
figure('units','normalized','outerposition',[0 0 1 1]);

%---------------------------------s2
subplot(2,2,4);hold all;
bar(z(1:3),coL2(1:3),'FaceColor','r');

bar(z(4:6),coL2(4:6),'FaceColor','b');
bar(z(7:9),coL2(7:9),'FaceColor','m');
bar(z(10:12),coL2(10:12),'FaceColor','g');
set(gca,'XTick',1:12);
set(gca,'XTickLabel',subs,'FontSize',8, 'FontWeight','bold');
grid;
title(' LQ R2','FontSize',12);
ylabel('% tot');
ylim([0 maL1]);

subplot(2,2,2);hold all;
bar(z(1:3),coH2(1:3),'FaceColor','r');

bar(z(4:6),coH2(4:6),'FaceColor','b');
bar(z(7:9),coH2(7:9),'FaceColor','m');
bar(z(10:12),coH2(10:12),'FaceColor','g');
set(gca,'XTick',1:12);
set(gca,'XTickLabel',subs,'FontSize',8, 'FontWeight','bold');

grid;
title(' HQ R2','FontSize',12);
ylabel('% tot');
ylim([0 maH1]);


%---------------------------------s1
subplot(2,2,3);hold all;
bar(z(1:3),coL1(1:3),'FaceColor','r');

bar(z(4:6),coL1(4:6),'FaceColor','b');
bar(z(7:9),coL1(7:9),'FaceColor','m');
bar(z(10:12),coL1(10:12),'FaceColor','g');
set(gca,'XTick',1:12);
set(gca,'XTickLabel',subs,'FontSize',8, 'FontWeight','bold');

grid;
title(' LQ R1', 'FontSize',12);
ylabel('% tot');
ylim([0 maL1]);

subplot(2,2,1);hold all;
bar(z(1:3),coH1(1:3),'FaceColor','r');

bar(z(4:6),coH1(4:6),'FaceColor','b');
bar(z(7:9),coH1(7:9),'FaceColor','m');
bar(z(10:12),coH1(10:12),'FaceColor','g');
set(gca,'XTick',1:12);
set(gca,'XTickLabel',subs,'FontSize',8, 'FontWeight','bold');

grid;
title(' HQ R1','FontSize',12);
ylabel('% tot');
%ylim([0 maa]);
ylim([0 maH1]);

end

%=================subfs
function [msh_1,msh_2, msl_1, msl_2,maL1,maH1]=metrics_per_cycle_plots_R1R2_HL_trim(msh_R1,msh_R2,msl_R1,msl_R2);

%--------------------analyses per cycle distribution of HQ substitutions
%          1. it should be flat across cycles
%             if it gradually grows in the end of read (as LQ) it should be reported 
%          2. it should be symmetrical between R1 and R2

%          it is trimmed (and reported) if last two cycles are stupidly inflated  
%           
%-------------------------OUTPUT
%msh_1,msh_2, msl_1, msl_2=

%  Normalised to all subs together (as percentage of all errors *100)
%-mim plots fr all subs, MSH MSL two samples R1 vs R2

%MSH=12 columns, ncyc lines (length n1 of each column)

%AC  AG AT CA CG CT GA GC GT TA TC TG
%1   2   3  4  5  6  7  8  9 10  11 12

%AC_dif=msh_2(:,1)-msh_1(:,1);
%TC_dif=msh_2(:,11)-msh_1(:,11);

ncyc1=length(msh_R1);
ncyc2=length(msh_R2);
n1=max(ncyc1,ncyc2);

%--------------------trim: zero nn1, nn2=2 so far
nn1=0;
nn2=2;
mshh_R1=msh_R1((nn1+1:ncyc1-nn2),:);
mshh_R2=msh_R2((nn1+1:ncyc2-nn2),:);

msll_R1=msl_R1((nn1+1:ncyc1-nn2),:);
msll_R2=msl_R2((nn1+1:ncyc2-nn2),:);

%--------------------continue and Normalise in a group HRi
all_R1=sum(sum(mshh_R1))+sum(sum(msll_R1));
all_R2=sum(sum(mshh_R2))+sum(sum(msll_R2));

all=all_R1+all_R2;

msh_1=100*mshh_R1/all;
msl_1=100*msll_R1/all;

msh_2=100*mshh_R2/all;
msl_2=100*msll_R2/all;


%1------------------------are there any ends of read effects to trim?
% ----------------pick up most important contributors (subs) and look at
% their tails (end of reads)
%sh1=sum(msh_1);
%ssh1=signif(sh1);

%-------------------if no slopes high positive (>2), no trimming
%[slopes_sub_vec,subs]= slopes_regressionSept21(msh_1);

%-------------------AC per cycle, HR1
%x=1:length(msh_1(:,1));
%[msh_1(:,1) 100*signif(msh_1(:,1))' x']




%============================PLOTS if needed

[maH1,maL1]=plots_MS(msh_1,msh_2, msl_1, msl_2);

end

%==============================subfunctions
function [maH1,maL1]=plots_MS(msh_1,msh_2, msl_1, msl_2)
%------------------------FOR PLOTS limits

     n1=length(msh_1(:,1));
     maH=5;
     maL=16;

maL1=max(max(max(msl_1),max(msl_2)));
%maL1=max(maL,mal);


maH1=max(max(max(msh_1),max(msh_2)));
%maH1=max(maH,mah);

n1=length(msh_1(:,1));% number of cycles

%=============================LQ

%figure;
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,4);
%plot(msl_2);
grid;
hold on;
plot(msl_2(:,1),'r');
plot(msl_2(:,2),'r');
plot(msl_2(:,3),'r');

plot(msl_2(:,4),'b');%C2A
plot(msl_2(:,5),'b');
plot(msl_2(:,6),'b');

plot(msl_2(:,7),'m');
plot(msl_2(:,8),'m');
plot(msl_2(:,9),'m');%G2T

plot(msl_2(:,10),'g');
plot(msl_2(:,11),'g');
plot(msl_2(:,12),'g');

plot(msl_2(:,1),'r-.','LineWidth',2);%AC
plot(msl_2(:,4),'b-.','LineWidth',2);%CA
plot(msl_2(:,9),'m-.','LineWidth',2);%AC
title('Low Q  R2');
xlabel('cycle number; AC(r) GT(m) CA(b)');
ylim([0 maL1]);
xlim([0 n1]);

subplot(2,2,3);
%plot(msl_1);
grid;
hold on;

plot(msl_1(:,1),'r');
plot(msl_1(:,2),'r');
plot(msl_1(:,3),'r');

plot(msl_1(:,4),'b');
plot(msl_1(:,5),'b');
plot(msl_1(:,6),'b');

plot(msl_1(:,7),'m');
plot(msl_1(:,8),'m');
plot(msl_1(:,9),'m');

plot(msl_1(:,10),'g');
plot(msl_1(:,11),'g');
plot(msl_1(:,12),'g');

plot(msl_1(:,1),'r-.','LineWidth',2);%AC
plot(msl_1(:,4),'b-.','LineWidth',2);%CA
plot(msl_1(:,9),'m-.','LineWidth',2);%AC
xlabel('cycle number; AC(r) GT(m) CA(b)');
title('Low Q  R1');
%xlabel('cycle number; AC(r) GT(m) CA(b)');
ylabel('%');
ylim([0 maL1]);
xlim([0 n1]);


%===============----------------------------HQ R1 R2-on top


subplot(2,2,2);
%plot(msh_2);
grid;
hold on;
plot(msh_2(:,1),'r');
plot(msh_2(:,2),'r');
plot(msh_2(:,3),'r');

plot(msh_2(:,4),'b');
plot(msh_2(:,5),'b');
plot(msh_2(:,6),'b');

plot(msh_2(:,7),'m');
plot(msh_2(:,8),'m');
plot(msh_2(:,9),'m');

plot(msh_2(:,10),'g');
plot(msh_2(:,11),'g');
plot(msh_2(:,12),'g');

plot(msh_2(:,1),'r-.','LineWidth',2);%AC
%plot(msh_2(:,2),'r','LineWidth',2);%AG
plot(msh_2(:,11),'g','LineWidth',2);%TC

%plot(msl_2(:,4),'b','LineWidth',2);%CA
%plot(msl_2(:,12),'g','LineWidth',2);%TG
plot(msh_2(:,4),'b-.','LineWidth',2);%CA
plot(msh_2(:,9),'m-.','LineWidth',2);%AC
xlabel('cycle number; AC(r) GT(m) CA(b)');

title('   HQ mismatch R2');
%xlabel('AG (r), TC(g), AC(r --)');
ylim([0 maH1]);
xlim([0 n1]);

subplot(2,2,1);
%plot(msh_1);
grid;
hold on;

plot(msh_1(:,1),'r');
plot(msh_1(:,2),'r');
plot(msh_1(:,3),'r');

plot(msh_1(:,4),'b');
plot(msh_1(:,5),'b');
plot(msh_1(:,6),'b');

plot(msh_1(:,7),'m');
plot(msh_1(:,8),'m');
plot(msh_1(:,9),'m');

plot(msh_1(:,10),'g');
plot(msh_1(:,11),'g');
plot(msh_1(:,12),'g');

plot(msh_1(:,1),'r-.','LineWidth',2);%AC
%plot(msh_1(:,9),'m','LineWidth',2);%GT
%plot(msh_1(:,2),'r','LineWidth',2);%AG
plot(msh_1(:,11),'g','LineWidth',2);%TC
plot(msh_1(:,4),'b-.','LineWidth',2);%CA
plot(msh_1(:,9),'m-.','LineWidth',2);%AC
xlabel('cycle number; AC(r) GT(m) CA(b)');


title('HQ mismatch plots for 12 subs, R1');
%xlabel('AG (r), TC(g), AC(r -.)');
ylabel('%');
ylim([0 maH1]);
xlim([0 n1]);
end

