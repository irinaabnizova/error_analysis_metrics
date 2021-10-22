function [predict_level]=prediction_C2A(thr_art,art_ox)
%5---------------predictions and prolly based on likelihood of art-ox
%-----------add hoc, and levels based on outputs of observations of C2A level by
%Jo and Kim (two different methods) which were strongly correlated with our
%GT =related metrics
%As far as we know, C2A evel values of Kim depends on sequence depth

%-------------INPUTS
%thr_art=0.8 : generally, high proportion of GT/CA to meanTi
%  art_ox(R1 and R2)  generally, proportion of GT/CA to meanTi
%-----------------------------OUTPUT
    %[predict_level]=predict_artefactOx(thr_art,art_ox,pred_low,pred_med,pred_high);
    % for R1 and R2
 
    
    predict_level=[];
    
    pred_high=3;
    pred_med=2;
    pred_low=1;
    pred_no=0;
    
    artR1=art_ox(1);
    artR2=art_ox(2);
    %was  somehow arbitrary (based on many examples) thr_art=0.8
     if artR1 >= thr_art | artR2 >= thr_art,
     %display('is likely C2A artefact');
     %art_ox
     [predict_level]=predict_artefactOx(thr_art,art_ox, pred_low,pred_med,pred_high);
 
     else
     %display('not likely C2A artefact');
     predict_level=[pred_no,pred_no];
     end
    
end
%===============================subfunction
function [predict]=predict_artefactOx(thr,art_ox, pred_low,pred_med,pred_high)
 %---------------computes predicted level of C2A based on its proportion to
 %nearest Ti (art_oxi) and observations from previous methods (should be
 %adjusted with time)
 
 %---------------------INPUTS
 %thr ~ 0.8: shows closeness enough to the Ti values (high contribution to
 %the error rate HQ substitutions)
 
 %art_ox = proportion of max(CA or GT) to nearest Ti (R1 and R2)
 %pred_low,pred_med,pred_high = observations from previous methods (should be
 %adjusted with time)
 
%---------------------initiate predictions for R1 R2
    predict=[0,0];
 
    %thr=0.8  ???
    thrM=thr+0.2;
    thrH=thr+0.6;
%---------------R1
    if art_ox(1)>=thrH,
         predict(1)=[pred_high];
    end
    if art_ox(1)>=thrM && art_ox(1)<thrH,
         predict(1)=[pred_med];
    end
    if art_ox(1)>=thr && art_ox(1)<thrM,
         predict(1)=[pred_low];
    end
 
 
%---------------R2
    if art_ox(2)>=thrH,
        predict(2)=[pred_high];
    end
    if art_ox(2)>=thrM && art_ox(2)<thrH,
        predict(2)=[pred_med];
    end
    if art_ox(2)>=thr && art_ox(2)<thrM,
        predict(2)=[pred_low];
    end
 
end