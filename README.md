# Generic Microsatellite Identification and Size Distribution

An analysis pipeline for identifying microsatellites from sequencing data. While designed for amplicon based sequencing projects, it could technically work for WES (though without necessary corrections for the higher false positive rate).

Pipeline requires:
1. Input paired FASTQs
    * gzipped
    * name formatted as *_R1.fastq.gz and *_R2.fastq.gz
2. BED file with loci of microsatellites of interest

The pipeline will produce histograms of microsatellite sizes as outputs for each sample.
