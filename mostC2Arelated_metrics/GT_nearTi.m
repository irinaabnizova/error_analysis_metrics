function [gt_nearTi]=GT_nearTi(ca1,gt1,ca2,gt2,tc1,ag1,tc2,ag2)
%--------GT/TC ratio to TC (the closest Ti) OR CA to AG
% ad hoc metric- just from frequent observations

%-----------------OUTPUTS
%  gt_nearTi=[max([gt_tc_1,ca_ag_1]),max([gt_tc_2,ca_ag_2])];

%-----------------INPUTS
%tc1,ag1,tc2,ag2- Ti pair nearest to GT/CA at the general profile plot
%ca1,gt1,ca2,gt2- CA and GT counts R1 and R2 

gt_nearTi=[];
%---------------------------body
    gt_tc_1=gt1/tc1;
    gt_tc_2=gt2/tc2;
    gt_tc=[gt_tc_1,gt_tc_2];
 
    ca_ag_1=ca1/ag1;
    ca_ag_2=ca2/ag2;
    ca_ag=[ca_ag_1,ca_ag_2];
        
    gt_nearTi=[max([gt_tc_1,ca_ag_1]),max([gt_tc_2,ca_ag_2])];
end