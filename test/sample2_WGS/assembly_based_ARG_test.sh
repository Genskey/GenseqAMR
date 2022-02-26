../../softs/ncbi-blast-2.9.0+/bin/blastn -evalue 1e-5 -num_threads 6 -outfmt 0 -num_alignments 10000 -query sample2.genome.fasta -db ../../DB/GenseqResDB/card.nucl -out sample2.AMR.m0
perl ../../bin/m0_format_convert.pl sample2.AMR.m0 sample2.AMR.m6
perl ../../bin/m6_filter2anno.pl sample2.AMR.m6 sample2.AMR.m6.filter.anno --long_query_seq --sub_coverage 0.4
perl ../../bin/ARG_detect2Anno.pl sample2.AMR.m6.filter.anno > sample2.AMR.detect.anno
perl ../../bin/obtain_candidate_feature.pl sample2.AMR.detect.anno > sample2.AMR.candidate_feature.stat
perl ../../bin/predict_AST.contig.pl sample2.AMR.candidate_feature.stat --pathogen "Klebsiella pneumoniae" > sample2.predict.xls
rm -rf sample2.AMR.m0 
