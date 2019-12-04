import "fastqc.wdl" as fastqc
import "single_sample_generic_microsatellite_identificaiton.wdl" as single_sample_generic_microsatellite_identificaiton

workflow GenericMicrosatelliteIdentification {
    # this workflow is a 'super' workflow

    # Input files
    Array[File] r1_files
    Array[File] r2_files

    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)

    String output_zip_name

    File loci_bed
    
    # Reference files
    File ref_fasta
    File ref_fasta_index
    File ref_index_1
    File ref_index_2
    File ref_index_3
    File ref_index_4
    File ref_rev_1
    File ref_rev_2

    # git repo info for report
    String git_repo_url
    String git_commit_hash

    scatter ( item in fastq_pairs ) {
        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = item.left
        }
        
        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = item.right
        }

        call single_sample_generic_microsatellite_identificaiton.SingleSampleGenericMicrosatelliteIdentification as single_sample_process {
            input:
                r1_fastq = item.left,
                r2_fastq = item.right,
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
    }



}