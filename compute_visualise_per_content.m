function [name_file_out1,name_file_out2]=compute_visualise_per_content(file_in,name)
         
       % computes, visualises and saves nucleotide content info
       
      % # nucleotide content_metrics each subs for tag, R1
      %substitution  12 13 14 21 23 24 31 32 34 41 42 43
      %significance  -0.90 0.40 -0.98 0.38 -0.82 1.97 0.69 -0.79 1.37 -0.71 0.29 -0.89
      %contribution  0 1 0 0 0 1 1 0 1 0 1 0
      %coef variation  1.62 1.82 1.25 1.87 1.27 1.04 1.31 1.28 1.58 1.86 2.00 1.49
      %relative range  0.45 0.52 0.42 0.69 0.39 0.26 0.37 0.37 0.53 0.53 0.55 0.41
      %max prev letter  2 2 2 2 2 2 2 2 4 2 2 3
      %max next letter  2 3 3 1 3 3 3 3 3 3 3 3
      %max window X.Y   22 23 23 41 11 23 23 23 43 23 23 33

       %--------------OUTPUTS: info about each of 12 susbs
       
       
       tags={'SET','RCH','RCL','PRCNH','PRCNL'};
       subss={'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG'};


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
         %[mshp_1,mshp_2, mslp_1, mslp_2,maL1,maH1]=metrics_per_cycle_plots_R1R2_HL_trim(msh_1,msh_2,msl_1,msl_2);
         %title([' Per cycle mismatch,  HQ R1, ',name1],'Interpreter','none','FontSize',14);
 
         
    %4-----PER NUCLEOTIDE CONTENT
    %4.1---------------------------------visualise
     [pr_subR1,n_subR1,winAX1,winCX1,winGX1,winTX1,ma1,mcw1,prAX1,prCX1,prGX1,prTX1]=metrics_nucleotides_pr_next_Ri(winHR1,E12(1));
          title([' Previous-Window-Next bases effect,  HQ R1, ',name1],'Interpreter','none','FontSize',14);

          %win_subR1=[winAX1;winCX1;winGX1;winTX1];% 12 columns of win X.Y for each subs
          
                 
     [pr_subR2,n_subR2,winAX2,winCX2,winGX2,winTX2,ma2,mcw2,prAX2,prCX2,prGX2,prTX2]=metrics_nucleotides_pr_next_Ri(winHR2,E12(2));
     title([' Previous-Window-Next bases effect,  HQ R2, ',name1],'Interpreter','none','FontSize',14);
     
     %4.2---------------------------------Ti content pu ch cv rel_range maxZ
          display('effect of previous and next nucleotide on Ti');
          ind_Ti=[2,6,7,11];
          [pr_eff_1i,pr_eff_2i,n_eff_1i,n_eff_2i]=effect_prev_next_nucl_noTi(ind_Ti,ind_Ti,pr_subR1,pr_subR2,n_subR1,n_subR2);
           
          display('signif prev next nucletide on Ti');
          thr_sig_cont=1.1
          [effect_prev_next_Ti]=effect_prev_next_nucleotide_Ti(thr_sig_cont,pr_subR1,pr_subR2,n_subR1,n_subR2);
       
        %4.2=======================all subs: effect window
        
        display('is there any systematic bias in windows X.Y frequency, all subs?');
          ind_R1=[1:12];
          ind_R2=[1:12];
   
          [win_eff_1,win_eff_2]=effect_window_max(ind_R1,ind_R2,winHR1,winHR2);
          [a]=plot_win_effect(prAX1,prCX1,prGX1,prTX1,prAX2,prCX2,prGX2,prTX2,win_eff_1,win_eff_2);

          [max_prev_1,max_prev_2]=max_prev_nucleotide(pr_subR1,pr_subR2);
          [max_next_1,max_next_2]=max_next_nucleotide(n_subR1,n_subR2);
          
          sub_pr_n_win_1=[max_prev_1;max_next_1(2,:);win_eff_1(4,:);win_eff_1(2,:);win_eff_1(3,:)];
          sub_pr_n_win_2=[max_prev_2;max_next_2(2,:);win_eff_2(4,:);win_eff_2(2,:);win_eff_2(3,:)];
          
          
          %-----------------------contributions of each subs in error rate
          thrS=1.1
          [sH12_con,indd1,indd2,sub_contrib1,sub_contrib2,con_12,sH12]=main_sub_contributors_ind(thrS,coH1,coH2);
          %sH12_con
          
          %------------------final info for subs
          
          sub_sig_con_cv_rra_pr_n_win_1=[max_prev_1(1,:);sH12_con(:,2)';con_12(:,2)';win_eff_1(2,:);win_eff_1(3,:);max_prev_1(2,:);max_next_1(2,:);win_eff_1(4,:)]
          sub_sig_con_cv_rra_pr_n_win_2=[max_prev_2(1,:);sH12_con(:,4)';con_12(:,3)';win_eff_2(2,:);win_eff_2(3,:);max_prev_2(2,:);max_next_2(2,:);win_eff_2(4,:)]
   

         %--------------------save main info about subs content
        [name_file_out1,name_file_out2]=save_content_tag(name,sub_sig_con_cv_rra_pr_n_win_1,sub_sig_con_cv_rra_pr_n_win_2);
        % saves sub_sig_con_cv_rra_pr_n_win_1 =
   
           %[name_out]=save_cycle_metrics_tag(name,infla_Ti_R12);
           %name_file_out= sprintf('%s_effect_nucleotide_metrics.txt',name);
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
function [pr1_subRi,n1_subRi,pr1A,pr1C,pr1G,pr1T,ma,mcw,prAX,prCX,prGX,prTX]=metrics_nucleotides_pr_next_Ri(new_winRi,toti);
%Normalised : by total Bases
%--------------------input :4    4x12  matrices: 1n1, 4n1
%new_winR1=winLR1m;
%tot1=totm;
%-----------------for 2 subs profiles: 4n1 1n1 (pr1 effect)
%z1=1:12;

%-last output  mcw=max_count_win;

%subs={'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG'};


%------------------ pr1 n1

[pr1_subRi,n1_subRi, group_nextRi] = win3_pr1_n1( new_winRi );


pe_pr1n=100*pr1_subRi/toti;
pe_n1n=100*n1_subRi/toti;

%pe2_pr1n=100*pr1_subR2/tot2;
%pe2_n1n=100*n1_subR2/tot2;

pr1A=new_winRi(1:4,:);% first four rows: prA just sum up
pr1C=new_winRi(5:8,:);
pr1G=new_winRi(9:12,:);
pr1T=new_winRi(13:16,:);


%------------------------plots

[pe_prA,pe_prC,pe_prG,pe_prT,pe_n1A,pe_n1C,pe_n1G,pe_n1T,prAX,prCX,prGX,prTX,ma,mcw]=plot_nucleotide_content(pe_pr1n,pe_n1n,new_winRi);


end

%=================sufun content plot
function [pe_prA,pe_prC,pe_prG,pe_prT,pe_n1A,pe_n1C,pe_n1G,pe_n1T,prAX,prCX,prGX,prTX,ma,mcw]=plot_nucleotide_content(pe_pr1n,pe_n1n,new_winRi)

%Normalised 15 Sept: by total Bases
%--------------------input :4    4x12  matrices: 1n1, 4n1
%new_winR1=winLR1m;
%tot1=totm;
%-----------------for 2 subs profiles: 4n1 1n1 (pr1 effect)
%z1=1:12;

%-last output  mcw=max_count_win;

subs={'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG'};


%display('===========DS SS and ALL SUBs 8371 H Pu');
%DS_SS_ALLF=[subNF' (subALL_F'-subNF') subALL_F']
 %[z]=plot_bar_DS_SS_ALL_num(DS_SS_ALLF);
 %set(gca,'XTick',1:12);set(gca,'XTickLabel',subs);
z=[0,0,0,0,0,0,0,0,0,0,0,0];
z4=[z' z' z' z'];% four columns, each length=12

%----------------------prs
        
   pe_prA = [ pe_pr1n(1,:)'  z'  z' z'];
    pe_prC = [ z' pe_pr1n(2,:)'  z' z'];
     pe_prG = [ z' z' pe_pr1n(3,:)'  z' ];
      pe_prT = [ z' z' z' pe_pr1n(4,:)' ];
   
      %---------------n1-s
      pe_n1A = [ pe_n1n(1,:)'  z'  z' z'];
    pe_n1C = [ z' pe_n1n(2,:)'  z' z'];
     pe_n1G = [ z' z' pe_n1n(3,:)'  z' ];
      pe_n1T = [ z' z' z' pe_n1n(4,:)' ];
       %--------------windows alphabetically
       %1. -------------pr1
pr1A=new_winRi(1:4,:);% first four rows: prA just sum up
pr1C=new_winRi(5:8,:);
pr1G=new_winRi(9:12,:);
pr1T=new_winRi(13:16,:);

max_count_win=max([max(pr1A),max(pr1C),max(pr1G),max(pr1T)]);
mcw=max_count_win;


prAX=[pr1A' z4 z4 z4];
prCX=[z4 pr1C' z4 z4];
prGX=[z4 z4 pr1G' z4];
prTX=[z4 z4 z4  pr1T' ];

    
  
      
      map=max(max([pe_pr1n]));
      man=max(max([pe_n1n]));
      ma=max(map, man);
%figure;
figure('units','normalized','outerposition',[0 0 1 1])


%---------------------------------s1


subplot(3,1,3);
    bar(pe_n1A,'r','EdgeColor','r');
    hold on;bar(pe_n1C,'b','EdgeColor','b');
    bar(pe_n1G,'m','EdgeColor','m');
    bar(pe_n1T,'g','EdgeColor','g');
%bar(pe_n1','grouped');
set(gca,'XTick',1:12);
%set(gca,'XTickLabel',subs,'FontSize',8);
set(gca,'XTickLabel',subs,'FontSize',10, 'FontWeight','bold');
grid;
title('Next base effect on substitutions','FontSize',12);
ylabel('% Next','FontSize',12);
%ylim([0 man]);
ylim([0 ma]);

subplot(3,1,2);
bar(prAX,'r','EdgeColor','r');
hold on;
bar(prCX,'b','EdgeColor','b');
bar(prGX,'m','EdgeColor','m');
bar(prTX,'g','EdgeColor','g');
%bar((new_winR1'));
set(gca,'XTick',1:12);
%set(gca,'XTickLabel',subs,'FontSize',8);
set(gca,'XTickLabel',subs,'FontSize',10, 'FontWeight','bold');
grid;
title(' Window effect on substitutions','FontSize',12);
ylabel('count win','FontSize',12);
%ylim([0 man]);


subplot(3,1,1);
    bar(pe_prA,'r','EdgeColor','r');
    hold on;bar(pe_prC,'b','EdgeColor','b');
    bar(pe_prG,'m','EdgeColor','m');
    bar(pe_prT,'g','EdgeColor','g');
set(gca,'XTick',1:12);
%set(gca,'XTickLabel',subs,'FontSize',8);
set(gca,'XTickLabel',subs,'FontSize',10, 'FontWeight','bold');
grid;
title(' Previous base effect on substitutions ','FontSize',12);
ylabel('% Prev','FontSize',12);
%ylim([0 map]);
ylim([0 ma]);

end

%=====================subfun
function [co_pr1,co_n1, group_next] = win3_pr1_n1( new_winR1 )

% computes one previous (pr) and one next (n) base from a 3-window;
% also computes reverse compliment window (16x12), reversing each column so it 
% can be summed up with its rev com subs, e.g. 

%    aCAa tGTt
%    aCAc gGCt
%    aCAg cGCt
%    aCAt aGCt

%-------------input 16 x12 : 16 windows X.Y for 12 subs(columns)
%winHR1m =

  %Columns 1 through 12
  %          AC            AG          AT            CA            CG   CT.
  %A.A    54702.00      59708.00      26088.00      39750.00      38233.00     104552.00      55283.00      37031.00      32088.00      27789.00
  %A.C    26519.00      42633.00      19370.00      23810.00      13718.00      41920.00      47190.00      15141.00      16428.00      23452.00
  %A.G    23918.00      55502.00      26800.00      14822.00      13979.00     107091.00      50619.00      16105.00      26936.00      33245.00
  %A.T    38847.00      89282.00      30491.00      25310.00      33900.00      79842.00      81134.00      37621.00      28114.00      30643.00
  %....
  
  %T.T    34606.00      60533.00      13077.00      33271.00      16103.00      51393.00      60683.00      17753.00       9813.00      10733.00

%------------------output
%co_n1=[nA;nC;nG;nT];%ACGT next  4 x 12 

%%co_pr1=[prA;prC;prG;prT];%ACGT next  4 x12

%%tgca for eact next T, G, C,A below: reverse complement with Illumina sequencing direction
%group_next=[reverse_row_matrix(n1T);reverse_row_matrix(n1G);reverse_row_matrix(n1C);reverse_row_matrix(n1A)];

    
  %------------------------------------windows
  subs={'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG'};
pairs_pr1={'AsA', 'AsC', 'AsG', 'AsT','CsA','CsC', 'CsG','CsT','GsA','GsC','GsG','GsT','TsA','TsC','TsG','TsT'};
pairs_n1={'AsT', 'CsT', 'GsT', 'TsT','AsG','CsG', 'GsG','TsG','AsC','CsC','GsG','TsC','AsA','CsA','GsA','TsA'};
subsP={'CA','CG','CT','TA','TC','TG'};
pairs_n1_rc={'TsT', 'GsT', 'CsT', 'AsT','TsG','GsG', 'CsG','AsG','TsC','GsC','CsG','AsC','TsA','GsA','CsA','AsA'};


%----------------------------computation
  %new_winR1=winLR2a;
  
  ma=max(max(new_winR1));
%1. -------------pr1
pr1A=new_winR1(1:4,:);% first four rows: prA just sum up
pr1C=new_winR1(5:8,:);
pr1G=new_winR1(9:12,:);
pr1T=new_winR1(13:16,:);

prA=sum(pr1A);
prC=sum(pr1C);
prG=sum(pr1G);
prT=sum(pr1T);

co_pr1=[prA;prC;prG;prT];

%coH_pr1_1a =

 % Columns 1 through 12

   %    8053.00      19904.00       9351.00       9259.00      10186.00      26187.00      17996.00
    %   7767.00      23936.00       7152.00       9196.00       4230.00      17585.00      23768.00
   %    4582.00      13320.00       6168.00       4709.00       4532.00      13943.00      12496.00
   %    6016.00      18774.00       7204.00       8758.00      10707.00      19380.00      20299.00

   
%coH_pr1_1a1 =

  %Columns 1 through 7

  %     8053.00      19904.00       9351.00       9259.00      10186.00      26187.00      17996.00
  %     7767.00      23936.00       7152.00       9196.00       4230.00      17585.00      23768.00
  %     4582.00      13320.00       6168.00       4709.00       4532.00      13943.00      12496.00
  %     6016.00      18774.00       7204.00       8758.00      10707.00      19380.00      20299.00


%2. ------------------next 1 n1
% next A 1,5,9,13: pr1=acgt
%next C  2,6,10,14
%next G 3,7,11,15
%T   4,8,12,16
n1A=[new_winR1(1,:);new_winR1(5,:);new_winR1(9,:);new_winR1(13,:)];% pr1 here are acgt
n1C=[new_winR1(2,:);new_winR1(6,:);new_winR1(10,:);new_winR1(14,:)];
n1G=[new_winR1(3,:);new_winR1(7,:);new_winR1(11,:);new_winR1(15,:)];
n1T=[new_winR1(4,:);new_winR1(8,:);new_winR1(12,:);new_winR1(16,:)];


nA=sum(n1A);
nC=sum(n1C);
nG=sum(n1G);
nT=sum(n1T);

co_n1=[nA;nC;nG;nT];%ACGT next
%group_next1=[n1T;n1G;n1C;n1A];%acgt vertically
%tgca below: reverse complement with Illumina sequencing direction
group_next=[reverse_row_matrix(n1T);reverse_row_matrix(n1G);reverse_row_matrix(n1C);reverse_row_matrix(n1A)];

end

%================subfun
function [newR_winR1]=reverse_row_matrix(new_winR1);

%input new_winR1  16 x 6 : 16 rows 6 columns

%--------------reverse all rows
    n=length(new_winR1(:,1));
    for i=1:n,
        k(i)=n-i+1;
        newR_winR1(i,:)=new_winR1(k(i),:);
         %newR_winR2(i,:)=new_winR2(k(i),:);
    end
end

%=============subfs on content Ti nonTi
function [pr_eff_1,pr_eff_2,n_eff_1,n_eff_2]=effect_prev_next_nucl_noTi(indd1,indd2,pr_subR1,pr_subR2,n_subR1,n_subR2)
      
%---------output
%effect_prev_next=[];
%%[ pr_eff;next_eff(2,:)]=
%effect_prev_next_Ti =

 %  Ti sub     13.00         24.00         31.00         42.00
 %  prev       2.00             0             0             0
 %  next          0             0          4.00          3.00


%------------------------INPUTs
%pr_subR2,n_subR2- counts of previous A C G T (lines) for each 12 subs
%(columns)

%pr_subR2 =

 % Columns 1 through 6
%         AC            AG            AT           CA ....
 % A    14976.00      47681.00      20618.00      37939.00      16316.00      62471.00
 % C    12501.00      70977.00      16512.00     302982.00      22012.00      55127.00
 % G    11336.00      37970.00      13966.00      51898.00      16036.00      47043.00
 % T    10229.00      41553.00      15210.00     297270.00      25260.00      51435.00

 % Columns 7 through 12

 %     45379.00      22684.00     201366.00      16299.00      67042.00      10420.00
 %     67488.00       8899.00      55241.00      15971.00      49281.00      11891.00
 %     34748.00      19026.00     164344.00      11389.00      29217.00       9475.00
  %    52550.00      29858.00     266549.00      18648.00      51054.00      14704.00


%n_subR2 = same format

%-------------------------------inndexes contributing subs: indd1 indd2
           
           %subs_num=subs_numbers;% 12 13  14 = AC AG AT etc
           
          % display('effect of previous and next nucleotide on contributing non Ti');
   
           %pr_subR1(:,indd1)=
%      37409.00     213746.00
%     297471.00      59005.00
%      48152.00     172785.00
%     290733.00     278015.00

% where indd1 =          4.00          9.00
                %3.3.2--------------other contributors pr next effect: just
          %output two max letters!
          
          [pr_eff_1,pr_eff_2]=effect_previous_letter(indd1,indd2,pr_subR1,pr_subR2);
          [n_eff_1,n_eff_2]=effect_next_letter(indd1,indd2,n_subR1,n_subR2);
          
end
%==========================subfunctions
function [n_eff_1,n_eff_2]=effect_next_letter(indd1,indd2,n_subR1,n_subR2)

%------------------------INPUTs

%----------indd1 indd2- per Read indexes in subs_num of inflated subs (not
%Ti)

%n_subR2,n_subR2- counts of previous A C G T (lines) for each 12 subs
%(columns)

%n_subR2 =

 % Columns 1 through 6
%         AC            AG            AT           CA ....
 % A    14976.00      47681.00      20618.00      37939.00      16316.00      62471.00
 % C    12501.00      70977.00      16512.00     302982.00      22012.00      55127.00
 % G    11336.00      37970.00      13966.00      51898.00      16036.00      47043.00
 % T    10229.00      41553.00      15210.00     297270.00      25260.00      51435.00

 % Columns 7 through 12

 %     45379.00      22684.00     201366.00      16299.00      67042.00      10420.00
 %     67488.00       8899.00      55241.00      15971.00      49281.00      11891.00
 %     34748.00      19026.00     164344.00      11389.00      29217.00       9475.00
  %    52550.00      29858.00     266549.00      18648.00      51054.00      14704.00


%n_subR2 = same format

%-------------------------------inndexes contributing subs: indd1 indd2
           
           subs_num=subs_numbers;% 12 13  14 = AC AG AT etc
           
           %display('effect of previous nucleotide on contributing non Ti');
   
           %n_subR1(:,indd1)=
%      37409.00     213746.00
%     297471.00      59005.00
%      48152.00     172785.00
%     290733.00     278015.00

         %3.3.2--------------other contributors n next effect: just
          %output two max letters!
          n_eff_1=[];
          if size(indd1)>0,
              %display('effect of next nucleotide on other contributing subs, R1');
              for k=1:length(indd1),
                 sub_type=subs_num(indd1(k));
                 %sig_prev_1=signif((n_subR1(:,indd1(:,k)))) 
                 [pu1,ch1,rang1,cv1,sig1]=pu_chas_rang_cv_sig_vec(n_subR1(:,indd1(:,k)));
                 mas1=max(sig1);
                 n_eff_1(:,k)=[subs_num(indd1(k)),pu1,ch1,rang1,cv1,mas1]';
              end
              
          end
          %n_eff_1
          
          n_eff_2=[];
          if size(indd2)>0,
              %('effect of previous and next nucleotide on other contributing subs, R2');
              for k=1:length(indd2),
                  sub_type=subs_num(indd2(k));
                  % sig_prev_2=signif((n_subR2(:,indd2(:,k))))
                 %[pu2,ch2,rang2,cv2]=purity_chastity_vec(n_subR2(:,indd2(:,k)));
                 [pu2,ch2,rang2,cv2,sig2]=pu_chas_rang_cv_sig_vec(n_subR2(:,indd2(:,k)));
                 mas2=max(sig2);
                 n_eff_2(:,k)=[subs_num(indd2(k)),pu2,ch2,rang2,cv2,mas2]';         
              end
          end
            %n_eff_2
end

%==============--------------sub2
function [pr_eff_1,pr_eff_2]=effect_previous_letter(indd1,indd2,pr_subR1,pr_subR2)

%------------------------INPUTs

%----------indd1 indd2- per Read indexes in subs_num of inflated subs (not
%Ti)

%pr_subR2,n_subR2- counts of previous A C G T (lines) for each 12 subs
%(columns)

%pr_subR2 =

 % Columns 1 through 6
%         AC            AG            AT           CA ....
 % A    14976.00      47681.00      20618.00      37939.00      16316.00      62471.00
 % C    12501.00      70977.00      16512.00     302982.00      22012.00      55127.00
 % G    11336.00      37970.00      13966.00      51898.00      16036.00      47043.00
 % T    10229.00      41553.00      15210.00     297270.00      25260.00      51435.00

 % Columns 7 through 12

 %     45379.00      22684.00     201366.00      16299.00      67042.00      10420.00
 %     67488.00       8899.00      55241.00      15971.00      49281.00      11891.00
 %     34748.00      19026.00     164344.00      11389.00      29217.00       9475.00
  %    52550.00      29858.00     266549.00      18648.00      51054.00      14704.00


%n_subR2 = same format

%-------------------------------inndexes contributing subs: indd1 indd2
           
           subs_num=subs_numbers;% 12 13  14 = AC AG AT etc
           
          % display('effect of previous nucleotide on contributing non Ti');
   
           %pr_subR1(:,indd1)=
%      37409.00     213746.00
%     297471.00      59005.00
%      48152.00     172785.00
%     290733.00     278015.00

         %3.3.2--------------other contributors pr next effect: just
          %output two max letters!
          pr_eff_1=[];
          if size(indd1)>0,
              %display('effect of previous nucleotide on other contributing subs, R1');
              for k=1:length(indd1),
                 sub_type=subs_num(indd1(k));
                 %sig_prev_1=signif((pr_subR1(:,indd1(:,k)))) 
                 [pu1,ch1,rang1,cv1,sig1]=pu_chas_rang_cv_sig_vec(pr_subR1(:,indd1(:,k)));
                 mas1=max(sig1);
                 pr_eff_1(:,k)=[subs_num(indd1(k)),pu1,ch1,rang1,cv1,mas1]';
              end
              
          end
         % pr_eff_1
          
          pr_eff_2=[];
          if size(indd2)>0,
              %display('effect of previous and next nucleotide on other contributing subs, R2');
              for k=1:length(indd2),
                  sub_type=subs_num(indd2(k));
                  % sig_prev_2=signif((pr_subR2(:,indd2(:,k))))
                 [pu2,ch2,rang2,cv2,sig2]=pu_chas_rang_cv_sig_vec(pr_subR2(:,indd2(:,k)));
                 mas2=max(sig2);
                 pr_eff_2(:,k)=[subs_num(indd2(k)),pu2,ch2,rang2,cv2,mas2]';         
              end
          end
           % pr_eff_2
end

function [pu,ch,rrang,cv,sig]=pu_chas_rang_cv_sig_vec(vec)
        %vec=layers(i,:);---positive values vector
        pu=-1;
        ch=-1;
        rrang=-1;
        cv=-1;  % unrealistic values if something is wrong
        %rrang=(amp-mi)/me;% relativenumber of means in a range
        %sig=signif(vec)
        sig=[];
        %-----------------------body
        mi=min(vec);
        me=mean(vec);
        sd=std(vec);
        if mi > 0,% all values are positive
        pu=max(vec)/sum(vec);
        amp=max(vec);
        ch=max(vec)/(max(vec)+second_max(vec));
        rrang=(amp-mi)/me;% number of means in a range
        cv=sd/me;
        sig=signif(vec);
        else
        display('negative values of data');
        %pu=-1;ch=-1;rang=-1;cv=-1;  
        end
end
%---------------------subfunction
function [second_max]=second_max(vecc)
   ma=max(vecc);
   vecc1=vecc;
   for i=1:length(vecc)
      if  vecc1(i)==ma,
        vecc1(i)=min(vecc);
      end
   end

   second_max=max(vecc1);
end

%======================
function [effect_prev_next_Ti]=effect_prev_next_nucleotide_Ti(thr_sig_cont,pr_subR1,pr_subR2,n_subR1,n_subR2)
      
%---------output
effect_prev_next_Ti=[];
%%[ pr_eff;next_eff(2,:)]=
%effect_prev_next_Ti =

 %  Ti sub     13.00         24.00         31.00         42.00
 %  prev       2.00             0             0             0
 %  next          0             0          4.00          3.00


%------------------------INPUTs
%pr_subR2,n_subR2- counts of previous A C G T (lines) for each 12 subs
%(columns)

%pr_subR2 =

 % Columns 1 through 6
%         AC            AG            AT           CA ....
 % A    14976.00      47681.00      20618.00      37939.00      16316.00      62471.00
 % C    12501.00      70977.00      16512.00     302982.00      22012.00      55127.00
 % G    11336.00      37970.00      13966.00      51898.00      16036.00      47043.00
 % T    10229.00      41553.00      15210.00     297270.00      25260.00      51435.00

 % Columns 7 through 12

 %     45379.00      22684.00     201366.00      16299.00      67042.00      10420.00
 %     67488.00       8899.00      55241.00      15971.00      49281.00      11891.00
 %     34748.00      19026.00     164344.00      11389.00      29217.00       9475.00
  %    52550.00      29858.00     266549.00      18648.00      51054.00      14704.00


%n_subR2 = same format

%-------------------------------inndexes for Ti subs
           ind_Ti=[2,6,7,11];
           subs_num=subs_numbers;% 12 13  14 = AC AG AT etc
           
           %display('effect of previous and next nucleotide on Ti');
   
           sig_prev_Ti_1=[signif(pr_subR1(:,2))' signif(pr_subR1(:,6))' signif(pr_subR1(:,7))' signif(pr_subR1(:,11))'];
           sig_next_Ti_1=[signif(n_subR1(:,2))' signif(n_subR1(:,6))' signif(n_subR1(:,7))' signif(n_subR1(:,11))'];
 
           sig_prev_Ti_2=[signif(pr_subR2(:,2))' signif(pr_subR2(:,6))' signif(pr_subR2(:,7))' signif(pr_subR2(:,11))'];
           sig_next_Ti_2=[signif(n_subR2(:,2))' signif(n_subR2(:,6))' signif(n_subR2(:,7))' signif(n_subR2(:,11))'];
          
           %thr_sig_cont=1.1
           
           %---------------------previous nucleotide effect
           pr_eff=[];
       for n=1:length(ind_Ti);% AG=subs_num(ind_Ti(1)) (ind_Ti=[2,6,7,11];
           vec=sig_prev_Ti_1(:,n);
           ns=0;eff_pr=0;
           for i1=1:length(vec),
               if vec(i1) >= thr_sig_cont,
                   ns=ns+1;
                   eff_pr=i1;
               end
           end
           pr_eff(:,n)=[subs_num(ind_Ti(n)),eff_pr]';
       end   
           %pr_eff
           
           %---------------------next nucleotide effect
           next_eff=[];
       for n=1:length(ind_Ti);% AG=subs_num(ind_Ti(1)) (ind_Ti=[2,6,7,11];
           vec=sig_next_Ti_1(:,n);
           ns=0;eff_next=0;
           for i1=1:length(vec),
               if vec(i1) >= thr_sig_cont,
                   ns=ns+1;
                   eff_next=i1;
               end
           end
           next_eff(:,n)=[subs_num(ind_Ti(n)),eff_next]';
       end   
           %next_eff
           
         effect_prev_next_Ti=[ pr_eff;next_eff(2,:)];
end

%=================subfun
function [z]=signif(y)

%-------------------------------computes Z score for a vector y

z=[];
me=mean(y);
sd=std(y);

    if sd>0,% there is variation
    
        for i=1:length(y),
            z(i)=(y(i)-me)/sd;
        end

    else
        z=[];
    end
    
end

%=========================effect window max
function [win_eff_1,win_eff_2]=effect_window_max(indd1,indd2,win_subR1,win_subR2)

%--------------------OUTPUTS
%win_eff_1(:,k)=[subs_num(indd1(k)),rang1,cv1,first_max_win]';
%------------------------INPUTs

%----------indd1 indd2- per Read indexes in subs_num of inflated subs (not
%Ti)

%pr_subR2,n_subR2- counts of previous A C G T (lines) for each 12 subs
%(columns)

%pr_subR2 =

 % Columns 1 through 6
%         AC            AG            AT           CA ....
 % A    14976.00      47681.00      20618.00      37939.00      16316.00      62471.00
 % C    12501.00      70977.00      16512.00     302982.00      22012.00      55127.00
 % G    11336.00      37970.00      13966.00      51898.00      16036.00      47043.00
 % T    10229.00      41553.00      15210.00     297270.00      25260.00      51435.00

 % Columns 7 through 12

 %     45379.00      22684.00     201366.00      16299.00      67042.00      10420.00
 %     67488.00       8899.00      55241.00      15971.00      49281.00      11891.00
 %     34748.00      19026.00     164344.00      11389.00      29217.00       9475.00
  %    52550.00      29858.00     266549.00      18648.00      51054.00      14704.00


%n_subR2 = same format

%-------------------------------inndexes contributing subs: indd1 indd2
           
           subs_num=subs_numbers;% 12 13  14 = AC AG AT etc
           wins=wins_numbers;
          % display('effect of previous nucleotide on contributing non Ti');
   
           %pr_subR1(:,indd1)=
%      37409.00     213746.00
%     297471.00      59005.00
%      48152.00     172785.00
%     290733.00     278015.00

         %3.3.2--------------other contributors pr next effect: just
          %output two max letters!
          win_eff_1=[];
          if size(indd1)>0,
              %display('effect of previous nucleotide on other contributing subs, R1');
              for k=1:length(indd1),
                 sub_type=subs_num(indd1(k));
                 %sig_prev_1=signif((pr_subR1(:,indd1(:,k)))) 
                 [pu1,ch1,rang1,cv1,sig1]=pu_chas_rang_cv_sig_vec(win_subR1(:,indd1(:,k)));
                 [first_max1,second_max1,indf1,inds1]=first_second_max_ind(win_subR1(:,indd1(:,k)));
                 %mas1=max(sig1);
                 win_eff_1(:,k)=[subs_num(indd1(k)),rang1,cv1,wins(indf1)]';
              end
              
          end
         % win_eff_1
          
          win_eff_2=[];
          if size(indd2)>0,
              %display('effect of previous and next nucleotide on other contributing subs, R2');
              for k=1:length(indd2),
                  sub_type=subs_num(indd2(k));
                  % sig_prev_2=signif((win_subR2(:,indd2(:,k))))
                 [pu2,ch2,rang2,cv2,sig2]=pu_chas_rang_cv_sig_vec(win_subR2(:,indd2(:,k)));
                 [first_max2,second_max2,indf2,inds2]=first_second_max_ind(win_subR1(:,indd1(:,k)));
                 %mas1=max(sig1);
                 win_eff_2(:,k)=[subs_num(indd2(k)),rang2,cv2,wins(indf2)]';
                 %win_eff_2(:,k)=[subs_num(indd2(k)),pu2,ch2,rang2,cv2,mas2]';         
              end
          end
           % win_eff_2
end
%====================first sec max
function [first_max,second_max,indf,inds]=first_second_max_ind(vecc)
 
%outputs vaues and indexes of first and second max

ma=max(vecc);
first_max=ma;
indf=0;
inds=0;
   vecc1=vecc;
   for i=1:length(vecc)
      if  vecc1(i)==ma,
        vecc1(i)=min(vecc);
        indf=i;
      end
   end

   second_max=max(vecc1);
       for i=1:length(vecc1)
          if  vecc1(i)==second_max,
          %vecc1(i)=min(vecc);
          inds=i;
          end
       end
  
end
%---------------------
function [a]=plot_win_effect(prAX1,prCX1,prGX1,prTX1,prAX2,prCX2,prGX2,prTX2,win_eff_1,win_eff_2)

 subss={'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG'};

a=1;
figure;
           subplot(2,1,1);
           plot(win_eff_1(3,:),'rp','LineWidth',2);
           title('cv(r) & relrange(m) win letter HR1');hold on;
           plot(win_eff_1(2,:),'mp','LineWidth',2);%title('relative range prev letter HR1');
  
           set(gca,'XTick',1:12);
           set(gca,'XTickLabel',subss,'FontSize',10, 'FontWeight','bold');
           grid;
           xlim([0 14]);
           subplot(2,1,2);
bar(prAX1,'r','EdgeColor','r');
hold on;
bar(prCX1,'b','EdgeColor','b');
bar(prGX1,'m','EdgeColor','m');
bar(prTX1,'g','EdgeColor','g');
set(gca,'XTick',1:12);grid;
set(gca,'XTickLabel',subss,'FontSize',10, 'FontWeight','bold');
%grid;
title(' Window effect on substitutions','FontSize',12);
ylabel('count win','FontSize',12);
xlim([0 14]);

figure;
           subplot(2,1,1);
           plot(win_eff_2(3,:),'bp','LineWidth',2);
           title('cv(b) & relrange(c) win letter HR2');hold on;
           plot(win_eff_2(2,:),'cp','LineWidth',2);%title('relative range prev letter HR1');
           set(gca,'XTick',1:12);
           set(gca,'XTickLabel',subss,'FontSize',10, 'FontWeight','bold');
           grid;
           xlim([0 14]);
           subplot(2,1,2);
bar(prAX2,'r','EdgeColor','r');
hold on;
bar(prCX2,'b','EdgeColor','b');
bar(prGX2,'m','EdgeColor','m');
bar(prTX2,'g','EdgeColor','g');grid;
set(gca,'XTick',1:12);
set(gca,'XTickLabel',subss,'FontSize',10, 'FontWeight','bold');
%grid;
title(' Window effect on substitutions','FontSize',12);
ylabel('count win','FontSize',12);
xlim([0 14]);      
          
end
%--------------------
function [lets]=letts_numbers

%unfold in vector of 4 letters A C G T
lets(1)=1;%a
lets(2)=2;% C
lets(3)=3;%G
lets(4)=4;%T
end

%=====================wins numbs
function [wins]=wins_numbers

%unfold in vector of 16 windows X.Y
wins(1)=11;%AA
wins(2)=12;% win_matr(2,1)=slopp(1);
wins(3)=13;%win_matr(3,1)=slopp(2);
wins(4)=14;%win_matr(4,1)=slopp(3);%

wins(5)=21;%win_matr(1,1)=slopp(4);
wins(6)=22;
wins(7)=23;%win_matr(3,2)=slopp(5);
wins(8)=24;%win_matr(4,2)=slopp(6);

wins(9)=31;%win_matr(1,2)=slopp(7);
wins(10)=32;%win_matr(2,2)=slopp(8);
wins(11)=33;
wins(12)=34;%win_matr(4,3)=slopp(9);

wins(13)=41;%win_matr(1,3)=slopp(10);
wins(14)=42;%win_matr(2,3)=slopp(11);
wins(15)=43;%win_matr(3,3)=slopp(12);
wins(16)=44;
end
%=========================
function [max_next_1,max_next_2]=max_next_nucleotide(n_subR1,n_subR2)
      % just get max prev and 
%---------output
max_next_1=[];
max_next_2=[];
%%[ n_eff;next_eff(2,:)]=
%max_next_next_Ti =

 %  Ti sub     13.00         24.00         31.00         42.00
 %  prev       2.00             0             0             0
 %  next          0             0          4.00          3.00


%------------------------INPUTs
%pr_subR2,n_subR2- counts of previous A C G T (lines) for each 12 subs
%(columns)

%pr_subR2 =

 % Columns 1 through 6
%         AC            AG            AT           CA ....
 % A    14976.00      47681.00      20618.00      37939.00      16316.00      62471.00
 % C    12501.00      70977.00      16512.00     302982.00      22012.00      55127.00
 % G    11336.00      37970.00      13966.00      51898.00      16036.00      47043.00
 % T    10229.00      41553.00      15210.00     297270.00      25260.00      51435.00

 % Columns 7 through 12

 %     45379.00      22684.00     201366.00      16299.00      67042.00      10420.00
 %     67488.00       8899.00      55241.00      15971.00      49281.00      11891.00
 %     34748.00      19026.00     164344.00      11389.00      29217.00       9475.00
  %    52550.00      29858.00     266549.00      18648.00      51054.00      14704.00


%n_subR2 = same format

%-------------------------------inndexes for Ti subs
           ind=[1:12];
           lets=letts_numbers;% 1 2 3 4
           subs_num=subs_numbers;
           
           display('max of previous nucleotide on subs');
           max_next_2=[];
           
              %display('effect of previous and next nucleotide on other contributing subs, R2');
              for k=1:length(ind),
                  sub_type=subs_num(ind(k));
                 [first_max2,second_max2,indf2,inds2]=first_second_max_ind(n_subR2(:,ind(:,k)));
                 max_next_2(:,k)=[subs_num(ind(k)),lets(indf2)]';
                        
              end
           
               max_next_1=[];
           
              %display('effect of previous and next nucleotide on other contributing subs, R2');
              for k=1:length(ind),
                  sub_type=subs_num(ind(k));
                 [first_max1,second_max1,indf1,inds1]=first_second_max_ind(n_subR1(:,ind(:,k)));
                 max_next_1(:,k)=[subs_num(ind(k)),lets(indf1)]';
                        
              end
end
function [max_prev_1,max_prev_2]=max_prev_nucleotide(pr_subR1,pr_subR2)
      % just get max prev and 
%---------output
max_prev_1=[];
max_prev_2=[];
%%[ pr_eff;next_eff(2,:)]=
%max_prev_next_Ti =

 %  Ti sub     13.00         24.00         31.00         42.00
 %  prev       2.00             0             0             0
 %  next          0             0          4.00          3.00


%------------------------INPUTs
%pr_subR2,n_subR2- counts of previous A C G T (lines) for each 12 subs
%(columns)

%pr_subR2 =

 % Columns 1 through 6
%         AC            AG            AT           CA ....
 % A    14976.00      47681.00      20618.00      37939.00      16316.00      62471.00
 % C    12501.00      70977.00      16512.00     302982.00      22012.00      55127.00
 % G    11336.00      37970.00      13966.00      51898.00      16036.00      47043.00
 % T    10229.00      41553.00      15210.00     297270.00      25260.00      51435.00

 % Columns 7 through 12

 %     45379.00      22684.00     201366.00      16299.00      67042.00      10420.00
 %     67488.00       8899.00      55241.00      15971.00      49281.00      11891.00
 %     34748.00      19026.00     164344.00      11389.00      29217.00       9475.00
  %    52550.00      29858.00     266549.00      18648.00      51054.00      14704.00


%n_subR2 = same format

%-------------------------------inndexes for Ti subs
           ind=[1:12];
           lets=letts_numbers;% 1 2 3 4
           subs_num=subs_numbers;
           
           display('max of previous nucleotide on subs');
           max_prev_2=[];
           
              %display('effect of previous and next nucleotide on other contributing subs, R2');
              for k=1:length(ind),
                  sub_type=subs_num(ind(k));
                 [first_max2,second_max2,indf2,inds2]=first_second_max_ind(pr_subR2(:,ind(:,k)));
                 max_prev_2(:,k)=[subs_num(ind(k)),lets(indf2)]';
                        
              end
           
               max_prev_1=[];
           
              %display('effect of previous and next nucleotide on other contributing subs, R2');
              for k=1:length(ind),
                  sub_type=subs_num(ind(k));
                 [first_max1,second_max1,indf1,inds1]=first_second_max_ind(pr_subR1(:,ind(:,k)));
                 max_prev_1(:,k)=[subs_num(ind(k)),lets(indf1)]';
                        
              end
end

%========================main subs contributions
function [sH12_con,indd1,indd2,sub_contrib1,sub_contrib2,cont_12,sH12]=main_sub_contributors_ind(thrS,coH1,coH2)
%----------------function main_contributors   
%-----------index them due to contribution: all Ti=1, Tv=0 only if
%comparable to Ti

indd1=[];indd2=[];sub_contrib1=[];sub_contrib2=[];
sH12=[];sH12_con=[];

[subs]=subs_numbers;

      %----------------------to detect main contributor subs, except Ti
           ind_Ti=[2,6,7,11];
           sH1=signif(coH1);
           sH2=signif(coH2);
           display('contribution of each substitution, R1 and R2');
           sH12=[sH1' sH2'];
       
           %initialise: all zero except Ti
           ind_con1=[0,1,0,0,0,1,1,0,0,0,1,0];
           ind_con2=[0,1,0,0,0,1,1,0,0,0,1,0]; 
           %--------------------get main HQ susbs contributor: per read
           %thrS=1.9;
           sub_contrib1=[];
           sub_contrib2=[];
           k1=0;k2=0;
           for i=1:length(sH1),
               if sH1(i) > thrS & i ~= 2 & i ~= 6 & i ~= 7 & i ~= 11,
                   k1=k1+1;
                   indd1(k1)=i;% indexes: which subs contribute
                   sub_contrib1(k1)=subs(i);
                   ind_con1(i)=1;
               end
              if sH2(i) > thrS & i ~= 2 & i ~= 6 & i ~= 7 & i ~= 11,
                   k2=k2+1;
                   sub_contrib2(k2)=subs(i);
                   indd2(k2)=i;% indexes: which subs contribute
                   ind_con2(i)=1;
              end
           end
           
        sH12_con=[subs' sH1' ind_con1' sH2' ind_con2'];
        cont_12=[subs' ind_con1' ind_con2'];
           
end
%==============
function [fn1,fn2]=save_content_tag(name,sub_sig_con_cv_rra_pr_n_win_1,sub_sig_con_cv_rra_pr_n_win_2)
%saves  Ti HQ cycle possible inflation metrics: 0 (No) or 1 (Yes), for
%Read1 (_R1) and Read2 (_R2)

%---------------------INPUT
% name = 38124/38124_2#154 (example, tag name with a path to directory name)

% sub_sig_con_cv_rra_pr_n_win_1 =
   
    % Columns 1 through 6
   
    %subs        12.00         13.00         14.00         21.00         23.00         24.00
    %signif      -0.90          0.40         -0.98          0.38         -0.82          1.97
    %contribut       0          1.00             0             0             0          1.00
    %cv           1.62          1.82          1.25          1.87          1.27          1.04
    %rel rang     0.45          0.52          0.42          0.69          0.39          0.26
    %maxPr        2.00          2.00          2.00          2.00          2.00          2.00
    %maxNex       2.00          3.00          3.00          1.00          3.00          3.00
    %max_win     22.00         23.00         23.00         41.00         11.00         23.00
   
   %  Columns 7 through 12
   
   %         31.00         32.00         34.00         41.00         42.00         43.00
   %          0.69         -0.79          1.37         -0.71          0.29         -0.89
   %          1.00             0          1.00             0          1.00             0
   %          1.31          1.28          1.58          1.86          2.00          1.49
   %          0.37          0.37          0.53          0.53          0.55          0.41
   %          2.00          2.00          4.00          2.00          2.00          3.00
   %          3.00          3.00          3.00          3.00          3.00          3.00
   %         23.00         23.00         43.00         23.00         23.00         33.00



%----------------OUTPUT
% name of an output file with a path: fn1 = sprintf('%s_.txt',name);

   ssccrpnw_1=sub_sig_con_cv_rra_pr_n_win_1;
   fn1 = sprintf('%s_content_metrics_R1.txt',name);
   dense_cr3=fopen(fn1,'w');
   
   fprintf(dense_cr3,'# nucleotide content_metrics each subs for tag, R1\n'); 
  
   fprintf(dense_cr3,'substitution  %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_1(1,1),ssccrpnw_1(1,2),ssccrpnw_1(1,3),ssccrpnw_1(1,4),ssccrpnw_1(1,5),ssccrpnw_1(1,6),ssccrpnw_1(1,7),ssccrpnw_1(1,8),ssccrpnw_1(1,9),ssccrpnw_1(1,10),ssccrpnw_1(1,11),ssccrpnw_1(1,12)); 
   fprintf(dense_cr3,'significance  %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',ssccrpnw_1(2,1),ssccrpnw_1(2,2),ssccrpnw_1(2,3),ssccrpnw_1(2,4),ssccrpnw_1(2,5),ssccrpnw_1(2,6),ssccrpnw_1(2,7),ssccrpnw_1(2,8),ssccrpnw_1(2,9),ssccrpnw_1(2,10),ssccrpnw_1(2,11),ssccrpnw_1(2,12)); 
   fprintf(dense_cr3,'contribution  %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_1(3,1),ssccrpnw_1(3,2),ssccrpnw_1(3,3),ssccrpnw_1(3,4),ssccrpnw_1(3,5),ssccrpnw_1(3,6),ssccrpnw_1(3,7),ssccrpnw_1(3,8),ssccrpnw_1(3,9),ssccrpnw_1(3,10),ssccrpnw_1(3,11),ssccrpnw_1(3,12)); 
   fprintf(dense_cr3,'coef variation  %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',ssccrpnw_1(4,1),ssccrpnw_1(4,2),ssccrpnw_1(4,3),ssccrpnw_1(4,4),ssccrpnw_1(4,5),ssccrpnw_1(4,6),ssccrpnw_1(4,7),ssccrpnw_1(4,8),ssccrpnw_1(4,9),ssccrpnw_1(4,10),ssccrpnw_1(4,11),ssccrpnw_1(4,12)); 
   fprintf(dense_cr3,'relative range  %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',ssccrpnw_1(5,1),ssccrpnw_1(5,2),ssccrpnw_1(5,3),ssccrpnw_1(5,4),ssccrpnw_1(5,5),ssccrpnw_1(5,6),ssccrpnw_1(5,7),ssccrpnw_1(5,8),ssccrpnw_1(5,9),ssccrpnw_1(5,10),ssccrpnw_1(5,11),ssccrpnw_1(5,12)); 
   fprintf(dense_cr3,'max prev letter  %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_1(6,1),ssccrpnw_1(6,2),ssccrpnw_1(6,3),ssccrpnw_1(6,4),ssccrpnw_1(6,5),ssccrpnw_1(6,6),ssccrpnw_1(6,7),ssccrpnw_1(6,8),ssccrpnw_1(6,9),ssccrpnw_1(6,10),ssccrpnw_1(6,11),ssccrpnw_1(6,12)); 
   fprintf(dense_cr3,'max next letter  %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_1(7,1),ssccrpnw_1(7,2),ssccrpnw_1(7,3),ssccrpnw_1(7,4),ssccrpnw_1(7,5),ssccrpnw_1(7,6),ssccrpnw_1(7,7),ssccrpnw_1(7,8),ssccrpnw_1(7,9),ssccrpnw_1(7,10),ssccrpnw_1(7,11),ssccrpnw_1(7,12)); 
   fprintf(dense_cr3,'max window X.Y   %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_1(8,1),ssccrpnw_1(8,2),ssccrpnw_1(8,3),ssccrpnw_1(8,4),ssccrpnw_1(8,5),ssccrpnw_1(8,6),ssccrpnw_1(8,7),ssccrpnw_1(8,8),ssccrpnw_1(8,9),ssccrpnw_1(8,10),ssccrpnw_1(8,11),ssccrpnw_1(8,12)); 

   fclose(dense_cr3);

 %F2===================

  ssccrpnw_2=sub_sig_con_cv_rra_pr_n_win_2;
   fn2 = sprintf('%s_content_metrics_R2.txt',name);
   dense_cr4=fopen(fn2,'w');
   
   fprintf(dense_cr4,'# nucleotide content_metrics each subs for tag, R2\n'); 
  
   fprintf(dense_cr4,'substitution  %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_2(1,1),ssccrpnw_2(1,2),ssccrpnw_2(1,3),ssccrpnw_2(1,4),ssccrpnw_2(1,5),ssccrpnw_2(1,6),ssccrpnw_2(1,7),ssccrpnw_2(1,8),ssccrpnw_2(1,9),ssccrpnw_2(1,10),ssccrpnw_2(1,11),ssccrpnw_2(1,12)); 
   fprintf(dense_cr4,'significance  %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',ssccrpnw_2(2,1),ssccrpnw_2(2,2),ssccrpnw_2(2,3),ssccrpnw_2(2,4),ssccrpnw_2(2,5),ssccrpnw_2(2,6),ssccrpnw_2(2,7),ssccrpnw_2(2,8),ssccrpnw_2(2,9),ssccrpnw_2(2,10),ssccrpnw_2(2,11),ssccrpnw_2(2,12)); 
   fprintf(dense_cr4,'contribution  %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_2(3,1),ssccrpnw_2(3,2),ssccrpnw_2(3,3),ssccrpnw_2(3,4),ssccrpnw_2(3,5),ssccrpnw_2(3,6),ssccrpnw_2(3,7),ssccrpnw_2(3,8),ssccrpnw_2(3,9),ssccrpnw_2(3,10),ssccrpnw_2(3,11),ssccrpnw_2(3,12)); 
   fprintf(dense_cr4,'coef variation  %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',ssccrpnw_2(4,1),ssccrpnw_2(4,2),ssccrpnw_2(4,3),ssccrpnw_2(4,4),ssccrpnw_2(4,5),ssccrpnw_2(4,6),ssccrpnw_2(4,7),ssccrpnw_2(4,8),ssccrpnw_2(4,9),ssccrpnw_2(4,10),ssccrpnw_2(4,11),ssccrpnw_2(4,12)); 
   fprintf(dense_cr4,'relative range  %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',ssccrpnw_2(5,1),ssccrpnw_2(5,2),ssccrpnw_2(5,3),ssccrpnw_2(5,4),ssccrpnw_2(5,5),ssccrpnw_2(5,6),ssccrpnw_2(5,7),ssccrpnw_2(5,8),ssccrpnw_2(5,9),ssccrpnw_2(5,10),ssccrpnw_2(5,11),ssccrpnw_2(5,12)); 
   fprintf(dense_cr4,'max prev letter  %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_2(6,1),ssccrpnw_2(6,2),ssccrpnw_2(6,3),ssccrpnw_2(6,4),ssccrpnw_2(6,5),ssccrpnw_2(6,6),ssccrpnw_2(6,7),ssccrpnw_2(6,8),ssccrpnw_2(6,9),ssccrpnw_2(6,10),ssccrpnw_2(6,11),ssccrpnw_2(6,12)); 
   fprintf(dense_cr4,'mex next letter  %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_2(7,1),ssccrpnw_2(7,2),ssccrpnw_2(7,3),ssccrpnw_2(7,4),ssccrpnw_2(7,5),ssccrpnw_2(7,6),ssccrpnw_2(7,7),ssccrpnw_2(7,8),ssccrpnw_2(7,9),ssccrpnw_2(7,10),ssccrpnw_2(7,11),ssccrpnw_2(7,12)); 
   fprintf(dense_cr4,'max window X.Y    %d %d %d %d %d %d %d %d %d %d %d %d\n',ssccrpnw_2(8,1),ssccrpnw_2(8,2),ssccrpnw_2(8,3),ssccrpnw_2(8,4),ssccrpnw_2(8,5),ssccrpnw_2(8,6),ssccrpnw_2(8,7),ssccrpnw_2(8,8),ssccrpnw_2(8,9),ssccrpnw_2(8,10),ssccrpnw_2(8,11),ssccrpnw_2(8,12)); 

   fclose(dense_cr4);
   
end