# error_analysis_metrics
metrics for raw substitution analysis after aligning Illumina reads: per cycle  and nucleotide content effect

![pic2git](https://user-images.githubusercontent.com/61786710/137915311-1e87ad8a-d444-46fc-a6fa-768c49fff630.png)

Therefore, to run it from a command line (i) for per cycle visualisation and metrics:

``name ='38912/38912_4#13';
thr_infla=1.3;
file_in = sprintf('%s_quality_error.txt',name);
[name_out]=compute_visualise_percycle(file_in,name,thr_infla)``

and (ii) for nucleotide content effect visualisation and metrics:

``name ='38912/38912_4#13'
file_in = sprintf('%s_quality_error.txt',name);
[name_file_outR1,name_file_outR2]=compute_visualise_per_content(file_in,name)``



*Required octave/Matlab*
