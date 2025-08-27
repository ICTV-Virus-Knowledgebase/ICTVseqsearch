#!/usr/bin/env bash
#
rsync --progress -hav \
	curtish@cheaha.rc.uab.edu:/data/project/ccts/curtish/ictv/VMR/VMR_to_BlastDB/blast \
	curtish@cheaha.rc.uab.edu:/data/project/ccts/curtish/ictv/VMR/VMR_to_BlastDB/blast_test \
	curtish@cheaha.rc.uab.edu:/data/project/ccts/curtish/ictv/VMR/VMR_to_BlastDB/processed_accessions_b.tsv \
	curtish@cheaha.rc.uab.edu:/data/project/ccts/curtish/ictv/VMR/VMR_to_BlastDB/processed_accessions_b.fa_names.tsv \
	.
