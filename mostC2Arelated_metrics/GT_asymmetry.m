function [symm_gt_ca]=GT_asymmetry(ca1,gt1,ca2,gt2)
%3.3-----------GT asymmetry
%here we look generally(in contrast to BROAD) at strand GT assymetry (sh be symmetry: #GT = #CA) 
    
%-------------------OUTPUT per R1 R2
%symm_gt_ca=[b1,b2];% sh be as close to zero as possible

    mi1=min([ca1,gt1]);
    mi2=min([ca2,gt2]);
    
    b1=-1;b2=-1;
    if mi1>0 & mi2>0,
    b1=(gt1-ca1)/mi1;% divided by smaller - not as in the BROAD paper
    b2=(gt2-ca2)/mi2;
    end
    
    symm_gt_ca=[b1,b2];% sh be as close to zero as possible
end