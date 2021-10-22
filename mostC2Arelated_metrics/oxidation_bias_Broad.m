function [bias_ox]=oxidation_bias_Broad(ca1,gt1,ca2,gt2)
%computes G-oxi as described by BROAD paper:
%%BROAD found CCG->CAG and GT more at R1, CA more at R2:
%they expect CA1<GT1 and CA2>GT2

%------------------output
%bias_ox=[bR1,bR2];

%-------------------------body of the function
    bR1=-1;bR2=-1;
    if ca1>0 & gt2>0,
    bR1=(gt1-ca1)/ca1;% divided by expected smaller- check the paper
    bR2=(gt2-ca2)/gt2;
    end
    bias_ox=[bR1,bR2];
end