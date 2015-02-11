## This is a Velvet de novo assembly pipeline
## Started 29 January 2015
## Needed for MLST software, ResFinder
## Needs to optimize for K, -exp_cov, and -cov_cutoff




for hashLen in 55
do

	for covCut in 50
	do

		for insLen in 0
		do

			for expCov in 40
			do

			echo "Hash Length:" $hashLen >> velvet_assemblies/summary.txt
			echo "Coverage Cutoff:" $covCut >> velvet_assemblies/summary.txt
			echo "Insert Length:" $insLen >> velvet_assemblies/summary.txt
			echo "Expected Coverage:" $expCov >> velvet_assemblies/summary.txt

			## MAKE velveth HASH TABLE
			velveth \
				velvet_assemblies \
				$hashLen \
				-fastq.gz \
				-shortPaired \
				-separate \
				/proj/julianog/sequence_reads/htsf_seq_backups/2014_11_24_Hajime_acineto_illumina/CLIN_A03_AGGCAGA-TAGATCG_L001_R1_001.fastq.gz \
				/proj/julianog/sequence_reads/htsf_seq_backups/2014_11_24_Hajime_acineto_illumina/CLIN_A03_AGGCAGA-TAGATCG_L001_R2_001.fastq.gz


			## DO THE ALIGNMENT PART
			velvetg \
				velvet_assemblies \
				-cov_cutoff $covCut \
				-ins_length $insLen \
				-exp_cov $expCov \
				| grep "Final graph has" >> velvet_assemblies/summary.txt

			done
		done
	done
done

