function [GT_meTi]=GT_averageTi(gt1,gt2,meTi)
%GT_meTi is similar to GT_Ti, but a little more general. 
%It is computed as a maximum between 
%(i) ratio of GT counts to mean(Ti) and 
%(ii) ratio CA to mean(Ti). -----------------not here yet!!!
%Ideally it is also less than 0.5.

%----also ad hoc: we observed that GT might be even more affected than CA

%-------------------INPUT
%gt1,gt2,meTi-2-vector

%----------------------OUTPUT
%GT_meTi

    GT_meTi=[];
  %3.5---------------------GT to meanTi
    GT_meTi1=gt1/meTi(1);
    GT_meTi2=gt2/meTi(2);
    GT_meTi=[GT_meTi1,GT_meTi2];
end