version 1.0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run FastViFi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
task RunFastViFi {
    input {
        File fastq1
        File? fastq2

        File Human_Reference
        File Virus_Reference
        File Kraken_db
      
        Int cpus
        Int preemptible
        String docker
        String sample_id
    }

    command <<<
        set -e


        export REFERENCE_REPO=`pwd`/viral_data
        export AA_DATA_REPO=`pwd`/data_repo
        export VIFI_DIR="/usr/local/src/ViFi"

        #~~~~~~~~~~~~~~~~~~~~~~~~
        # Untar the references  
        #~~~~~~~~~~~~~~~~~~~~~~~~
        tar -xvf ~{Human_Reference}
        tar -xvf ~{Virus_Reference}
        tar -xvf ~{Kraken_db}

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
         fastq1=~{fastq1}
         fastq2=~{fastq2}
      fi
      

      #~~~~~~~~~~~~~~~~~~~~~~~
      # run FastViFi
      #~~~~~~~~~~~~~~~~~~~~~~~

      python /usr/local/src/FastViFi/run_kraken_vifi_pipeline.py \
        --output-dir `pwd`/outdir \
        --input-file $fastq1 \
        --input-file-2 $fastq2 \
        --level sample-level-validation-intermediate \
        --kraken-path /usr/local/bin/kraken2 \
        --kraken-db-path `pwd`/kraken_datasets \
        --vifi-path /usr/local/src/ViFi/scripts/run_vifi.py \
        --virus hpv --human-chr-list /usr/local/src/FastViFi/test/human_chr_list.txt \
        --skip-bwa-filter \
        --keep-intermediate-files \
        --vifi-human-ref-dir /data/data_repo \
        --vifi-viral-ref-dir /data/viral_data \
        --docker
      

      mv outdir/output_hpv.clusters.txt ~{sample_id}.fastvifi.output_hpv.clusters.txt
      mv outdir/output_hpv.clusters.txt.range ~{sample_id}.fastvifi.output_hpv.clusters.txt.range
      
    >>>

    output {
        File output_clusters_txt_range="~{sample_id}.fastvifi.output_hpv.clusters.txt.range"
        File output_clusters_txt = "~{sample_id}.fastvifi.output_hpv.clusters.txt"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(Virus_Reference, "GB") + size(Human_Reference, "GB") +  size(Kraken_db, "GB") +  size(fastq1, "GB")*6 + 100) + " HDD"
        docker: docker
        cpu: cpus
        memory: "100GB"
    }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

workflow FastViFi {
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
        # Directories 
        #~~~~~~~~~~~~
        File Virus_Reference
        File Human_Reference
        File Kraken_db
      
        #~~~~~~~~~~~~
        # general runtime settings
        #~~~~~~~~~~~~
        Int preemptible = 2
        String docker = "brownmp/vifi:devel"

        

    }

    parameter_meta {
        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        cpus:{help:"CPU count"}
        docker:{help:"Docker image"}
    }


    #########################
    # run using given references 
    #########################
    call RunFastViFi{
        input:
            fastq1 = left,
            fastq2 = right,

            Human_Reference = Human_Reference,
            Virus_Reference = Virus_Reference,
            
            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id
    }
}


