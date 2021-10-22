function [art_ox,artR1,artR2]=artifact_oxi(ca1,gt1,ca2,gt2,meTi)
%computes==art_ox=proportion of max(gt or ca) to mean Ti within Ri
    
%--------------------OUTPUT
%art_ox=[artR1,artR2];- 2-vector
%---------------------INPUTS
%ca1,gt1,ca2,gt2
%meTi-average ti for R1 and R2

%----------------------body

    art_ox=[];

    ma1=max([ca1,gt1]);
    ma2=max([ca2,gt2]);
 
    artR1=ma1/meTi(1);% prop max(gt,ca) in Ti: the higher the more likely artC2A
    artR2=ma2/meTi(2);
    art_ox=[artR1,artR2];
    
