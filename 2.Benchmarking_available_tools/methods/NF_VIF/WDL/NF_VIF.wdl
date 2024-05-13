version 1.0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run NF-VIF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

task Run_NF_VIF {
    input {
        File fastq1
        File? fastq2

        File HPV_BWT2_index_tar
        File Human_Reference_tar
        File GRCh38_genome_2bit
        File ref_genome_fa
        File HPVs_db_fasta
        
        Int cpus
        Int preemptible
        String docker
        String sample_id
      }


    command <<<
        
       set -ex

       #~~~~~~~~~~~~~~~~~~~~~~~~
        # Untar the references  
        #~~~~~~~~~~~~~~~~~~~~~~~~
        tar -xvf ~{Human_Reference_tar} -C .
        rm -f ~{Human_Reference_tar}  
        tar -xvf ~{HPV_BWT2_index_tar} -C .
        rm -f ~{HPV_BWT2_index_tar}
       
        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]]
        then
            mkdir fastq
            tar -xvf ~{fastq1} -C fastq
            rm ~{fastq1}
            #fastqs=$(find fastq -type f)
            fastqs=($(pwd)/fastq/*)
            fastq1="${fastqs[0]}"
            fastq2="${fastqs[1]}"
        else 
            fastq1="~{fastq1}"
            fastq2="~{fastq2}" 
        fi

      # nf-vif is very particular about filenaming
      mkdir input_reads
      if [ "${fastq1: -3}" == ".gz" ]; then
        mv $fastq1 input_reads/reads_R1.fastq.gz
        mv $fastq2 input_reads/reads_R2.fastq.gz
      else
        mv $fastq1 input_reads/reads_R1.fastq
        mv $fastq2 input_reads/reads_R2.fastq
      fi
      
      
       nextflow /usr/local/src/nf-VIF/main.nf \
         --reads  "`pwd`/input_reads/*{1,2}.fastq*" \
         --genome hg38 \
         --bwt2_index `pwd`/human_reference \
         --fasta ~{ref_genome_fa} \
         --blatdb ~{GRCh38_genome_2bit} \
         --fasta_hpv ~{HPVs_db_fasta} \
         --bwt2_index_hpv  `pwd`/HPV_BWT2_index \
         --outdir nf_vif \
         --split_report --skip_trimming --skip_fastqc --skip_multiqc --nb_geno 1

      mkdir ~{sample_id}.nf-vif
      
      mv nf_vif/hpv_mapping/blat/*table*.csv ~{sample_id}.nf-vif/
      tar -zcvf ~{sample_id}.nf-vif.tar.gz ~{sample_id}.nf-vif/

       

      >>>


       output {
        File output_file="~{sample_id}.nf-vif.tar.gz"
      }


       runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(3*size(Human_Reference_tar, "GB") + size(fastq1, "GB")*6 + 100) + " HDD"
        docker: docker
        cpu: cpus
        memory: "100GB"
    }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

workflow NF_VIF {
    input {


        #~~~~~~~~~~~~
        # Sample ID
        #~~~~~~~~~~~~
        String sample_id
      
        #~~~~~~~~~~~~
        # FASTQ Files
        #~~~~~~~~~~~~
        File left
        File? right

        #~~~~~~~~~~~~
        # CPU count 
        #~~~~~~~~~~~~
        Int cpus = 10

        #~~~~~~~~~~~~
        # Resources
        #~~~~~~~~~~~~
        File HPV_BWT2_index_tar
        File Human_Reference_tar
        File GRCh38_genome_2bit
        File ref_genome_fa
        File HPVs_db_fasta
      
        #~~~~~~~~~~~~
        # general runtime settings
        #~~~~~~~~~~~~
        Int preemptible = 1
        String docker = "trinityctat/nf_vif:devel"

      
        

    }
 parameter_meta {
        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        docker:{help:"Docker image"}
        
    }


    #########################
    # run using given references 
    #########################
    call Run_NF_VIF{
        input:
            fastq1 = left,
            fastq2 = right,

            HPV_BWT2_index_tar = HPV_BWT2_index_tar,
            Human_Reference_tar = Human_Reference_tar,
            GRCh38_genome_2bit = GRCh38_genome_2bit,
            ref_genome_fa = ref_genome_fa,
            HPVs_db_fasta = HPVs_db_fasta,
      
            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id
    }
}

