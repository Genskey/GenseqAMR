perl ../../bin/convert_fq2fa.pl sample1.patho.fq.gz sample1.fa 1
../../softs/ncbi-blast-2.9.0+/bin/blastn -query sample1.fa -db ../../DB/GenseqResDB/card.nucl -evalue 1e-5  -num_threads 2 -outfmt 0 -num_alignments 10000  -out sample1.m0
perl ../../bin/m0_format_convert.pl sample1.m0 sample1.m6
perl ../../bin/m6_filter2anno.pl sample1.m6 sample1.m6.filter.anno
perl ../../bin/ARGdetect_stat.pl sample1.m6.filter.anno sample1.stat.xls sample1.reads.anno.lca sample1.reads.anno.lca.exact 2> arg_stat.log
perl ../../bin/ARGdetect_stat.filter.pl sample1.stat.xls > sample1.stat.filter.xls
#rm -rf sample1.m0

../../softs/minimap2-2.17_x64-linux/minimap2 -x sr -a -t 8 --secondary=no -L ../../DB/target_pathogens_Refgenome/Ref_genome.fasta sample1.fa > sample1.sam
perl ../../bin/calculate_cover_depth.pl sample1.sam --prefix sample1
#rm -rf sample1.sam

perl ../../bin/patho_ARGdetect.pl sample1.stat.filter.xls sample1.cvgstat > sample1.DrugClass_detectARGs.xls 2> arg_locate.log
perl ../../bin/predict_AST.read.pl --mNGSread sample1.DrugClass_detectARGs.xls >  sample1.predict.xls

