# T2T_IntSVTR
This repository contains the primary analysis scripts for the comparative evaluation of variant datasets from second- and third-generation sequencing technologies.

## scripts
以下是script目录下的脚本说明

| script | Description |
|-----------|-------------|
| svcmp.1kinv.py | Script for comparing INV variants in second- and third-generation datasets |
| svcmp.1kvntr_nwtr.py | Script for comparing TR variants in second- and third-generation datasets |
| svcmp.py | Script for comparing SV variants in second- and third-generation datasets |
| get_1019sv_nwsv_nwme_allsvtype_detail_pic_over05_svlenge50.rm05mesv.rm05trsv_splitdup.py | Scripts for calculating the total count of various types of SVs and for performing concordance analysis. |
| get_1019sv_nwsv_nwme_allsvtype_region_ratio_detail_pic_over05_svlenge50.rm05mesv.rm05trsv_splitdup_new_addratio.py | Script for analyzing the consistency of various SVs across different genomic regions |
| get_65allsv_1019sv_nwsv_nwme_allsvtype_detail_pic_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.precision.line.py | Scripts for precision analysis of various SVs |
| get_pic_1kallsv_nyuwaallsv_65_over05_svlen_distribution_svlenge50.rm05mesv.rm05trsv_splitdup_precision.line.py | Scripts for analyzing the precision of various types of SVs across different length ranges |
| get_65allsv_1019sv_nwsv_nwme_allsvtype_region_detail_pic_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.precision.py | Scripts for analyzing the precision of various types of SVs across different genomic regions |
| get_65allsv_1019sv_nwsv_nwme_allsvtype_detail_pic_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.recall.py | Scripts for recall analysis of various SVs |
| get_pic_1kallsv_nyuwaallsv_65_over05_svlen_distribution_svlenge50.rm05mesv.rm05trsv_splitdup_recall.py | Scripts for analyzing the recall of various types of SVs across different length ranges |
| get_65allsv_1019sv_nwsv_nwme_allsvtype_region_detail_pic_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.recall.py | Scripts for analyzing the recall of various types of SVs across different genomic regions |
| get_1019svtr_1019vntr_nwtr_reflen_region_ratio_detail_pic_over05_adjust_new.py | Scripts for analyzing the consistency of TRs associated with various reflens across different genomic regions |
| get_1019svtr_1019vntr_nwtr_motif_region_ratio_detail_pic_over05_adjust_new.py | Scripts for analyzing the consistency of TRs associated with various motifs across different genomic regions |
| get_pic_1ktr_nwtr_motif_recall_base05_or05.macmotif.rmnoreflen.py | Scripts for analyzing the recall of various motifs within TRs across different genomic regions |
| get_pic_1ktr_nwtr_reflen_recall_base05_or05.macmotif.rmnoreflen.py | Scripts for analyzing the recall of various reflens within TRs across different genomic regions |
| get_pic_1ktr_nwtr_reflen_precision_base05_or05.macmotif.rmnoreflen.rmonedata.py | Scripts for analyzing the precision of various reflens within TRs across different genomic regions |
| get_pic_1ktr_nwtr_motif_precision_base05_or05.macmotif.rmnoreflen.rmonedata.py | Scripts for analyzing the precision of various motifs within TRs across different genomic regions |
