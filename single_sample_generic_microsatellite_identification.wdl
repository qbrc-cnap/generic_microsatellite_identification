import "bwa_align.wdl" as bwa_align

workflow SingleSampleGenericMicrosatelliteIdentification {
    File r1_fastq
    File r2_fastq
    String basename = basename(r1_fastq, "_R1.fastq.gz")

    File loci_bed
    
    File ref_fasta
    File ref_fasta_index
    File ref_index_1
    File ref_index_2
    File ref_index_3
    File ref_index_4
    File ref_rev_1
    File ref_rev_2

    call remove_adapters {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            basename = basename
    }

    call merge_fastq_pairs {
        input:
            r1_fastq = remove_adapters.r1_fastq,
            r2_fastq = remove_adapters.r2_fastq,
            basename = basename
    }

    call convert_FQ_to_FA {
        input:
            fastq = merge_fastq_pairs.merged_fastq,
            basename = basename
    }

    call identify_with_phobos {
        input:
            fasta = convert_FQ_to_FA.fasta,
            basename = basename,
            loci_bed = loci_bed,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_index_1 = ref_index_1,
            ref_index_2 = ref_index_2,
            ref_index_3 = ref_index_3,
            ref_index_4 = ref_index_4,
            ref_rev_1 = ref_rev_1,
            ref_rev_2 = ref_rev_2
    }

    call bwa_align.perform_align as alignment {
        input:
            r1_fastq = .r1_fastq,
            r2_fastq = .r2_fastq,
            sample_name = basename,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            ref_bwt = ref_bwt,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_pac = ref_pac,
            ref_sa = ref_sa
    }
}

task remove_adapters {
    File r1_fastq
    File r2_fastq
    String basename

    command {
        cutadapt \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -o ${basename}_1.fastq \
            -p ${basename}_2.fastq \
            ${r1_fastq} \
            ${r2_fastq} > ${basename}.cutadapt_stats.txt;
    }

    output {
        File r1_fastq = "${basename}_1.fastq"
        File r2_fastq = "${basename}_2.fastq"
        File cutadapt_metrics = "${basename}.cutadapt_stats.txt"
    }
}

task merge_fastq_pairs {
    File r1_fastq
    File r2_fastq
    String basename

    command {
        /opt/software/bbmap/bbmerge.sh \
            in=${r1_fastq} \
            in2=${r2_fastq} \
            out=${basename}.fastq.gz \
            outu1=${basename}_1.unmerged.fastq.gz \
            outu2=${basename}_2.unmerged.fastq.gz \
        2> ${basename}.merge_stats.txt;
    }

    output {
        File merged_fastq = "${basename}.fastq.gz"
        File r1_unmerged_fastq = "${basename}_1.unmerged.fastq.gz"
        File r2_unmerged_fastq = "${basename}_2.unmerged.fastq.gz"
        File bbmerge_metrics = "${basename}.merge_stats.txt"
    }

}

task convert_fq_to_fa {
    File fastq
    String basename

    command {
        /opt/software/seqtk-1.3/seqtk seq \
            -a ${fastq} \
        > ${basename}.fasta
    }

    output {
        File fasta = "${basename}.fasta"
    }
}

task identify_with_phobos {
    File fasta
    String basename
    File loci_bed
    File ref_fasta
    File ref_fasta_index
    File ref_index_1
    File ref_index_2
    File ref_index_3
    File ref_index_4
    File ref_rev_1
    File ref_rev_2

    command {
        # Identify microsatellites with PHOBOS
        # flank 
        # --minScore 5 - above default of 0
        # --minLength_b 6 - 
        # --minUnitLen 1 - minimum repeating units
        # --maxUnitLen 6 - maximum repeating units
        # --flanking - print 150bp flanking nucleotides in read
        /opt/software/phobos-v3.3.12-linux/bin/phobos_64_libstdc++6 \
            ${fasta} \
            ${basename}.phobos \
            -M imperfect \
            --minScore 5 \
            --minLength_b 6 \
            --minUnitLen 1 \
            --maxUnitLen 6 \
            --flanking 150 \
            --outputFormat 3;
        # Identify largest microsatellite in read from PHOBOS output
        # Then outputs the data as interleaved FASTA to stdout.
        # The output is piped to shell scripts to separate the FASTA into
        # paired FASTA files
        python3 /opt/software/parse_phobos.py \
            --minflanking 10 \
            $f \
        | paste - - - - - - - - \
        | tee >(cut -f 1-2 | tr "\t" "\n" > ${basename}.ms_flanking.1.fa) \
        | cut -f 3-4 | tr "\t" "\n" > ${basename}.ms_flanking.2.fa;
        # Align the flanking reads to genome
        # Filter for good map quality, locus of interest, and concordant mapping
        # Collate read names
        /opt/software/bowtie2-2.3.5.1-linux-x86_64/bowtie2 \
            --no-unal \
            -f \
            -x ${ref_fasta} \
            -1 ${basename}.ms_flanking.1.fa \
            -2 ${basename}.ms_flanking.2.fa \
        2> ${basename}.aln_stats.txt \
        | /opt/software/samtools/bin/samtools view -q 10 -f 3 -L ${loci_bed} \
        | cut -f1 \
        | sed 's/aaa/\|/g' \
        | cut -f1 -d'|' \
        | sort \
        | uniq \
        > ${basename}.flt.names.txt;
        # Filter the PHOBOS output for reads in the filtered 'good' selection
        # Filter the FASTQ with seqtk, then run PHOBOS again
        /opt/software/seqtk-1.3/seqtk subseq \
            ${basename}.fasta ${basename}.flt.names.txt \
        > ${basename}.flt.fasta;
        /opt/software/phobos-v3.3.12-linux/bin/phobos_64_libstdc++6 \
            ${basename}.flt.fasta  \
            ${basename}.flt.phobos \
            -M imperfect \
            --minScore 5 \
            --minLength_b 6 \
            --minUnitLen 1 \
            --maxUnitLen 6 \
            --flanking 150 \
            --outputFormat 3;
        # Filter the PHOBOS with the python script again (with fltonly flag)
        python3 /opt/software/parse_phobos.py \
            --minflanking 10 \
            --fltonly \
            $f \
        | tail -n +30 \
        | awk -F '\t' '{print $6}' | sort | uniq -c | sort -nr \
        > ${basename}.ms_counts.txt
    }

    output {
        File microsatellite_counts = "${basename}.ms_counts.txt"
    }
}
