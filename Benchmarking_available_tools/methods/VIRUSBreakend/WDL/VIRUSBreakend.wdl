version 1.0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run VIRUSBreakend
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
task RunVIRUSBreakend {
    input {
        File fastq1
        File? fastq2

        # References
        File Human_Reference
        File Virus_Reference
        

        Int cpus
        Int preemptible
        String docker
        String sample_id
        Int disk
    }

    command <<<
        set -e



        #~~~~~~~~~~~~~~~~~~~~~~~~
        # Untar the references  
        #~~~~~~~~~~~~~~~~~~~~~~~~
        tar -xvf ~{Human_Reference}
        tar -xvf ~{Virus_Reference}
        


        #~~~~~~~~~~~~~~~~~~~~~~~~
        # Run VIRUSBreakend
        #   two specific cases 
        #       with two fastqs tared togeter and seperate 
        #~~~~~~~~~~~~~~~~~~~~~~~~
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


            bwa mem BWA_index/ref_genome.fa \
                -t ~{cpus} \
                $fastq1 \
                $fastq2 \
                > bwa_align.bam

            samtools view -bS bwa_align.sam > bwa_align.bam
            samtools sort bwa_align.bam -o bwa_sorted.bam

            #~~~~~~~~~~~~~~~~~~~~~~~
            # Run VIRUSBreakend
            #~~~~~~~~~~~~~~~~~~~~~~~
            virusbreakend \
                --kraken2db virusbreakenddb_20210401 \
                --output ~{sample_id}.virusbreakend.vcf \
                --reference `pwd`/BWA_index/ref_genome.fa \
                bwa_sorted.bam  \
                2>&1 | tee output_log_subset.txt
           
        
        else 
            #~~~~~~~~~~~~~~~
            # Run Alignment
            #~~~~~~~~~~~~~~~
            
            bwa mem BWA_index/ref_genome.fa \
                -t ~{cpus} \
                ~{fastq1} \
                ~{fastq2} \
                > bwa_align.sam

            samtools view -bS bwa_align.sam > bwa_align.bam
            samtools sort bwa_align.bam -o bwa_sorted.bam

            #~~~~~~~~~~~~~~~~~~~~~~~
            # Run VIRUSBreakend
            #~~~~~~~~~~~~~~~~~~~~~~~
            virusbreakend \
                --kraken2db virusbreakenddb_20210401 \
                --output ~{sample_id}.virusbreakend.vcf \
                --reference `pwd`/BWA_index/ref_genome.fa \
                bwa_sorted.bam  \
                2>&1 | tee output_log_subset.txt
        fi

        #~~~~~~~~~~~~~~~~~~~~~~~~
        # Tar the output
        #~~~~~~~~~~~~~~~~~~~~~~~~
        # tar -czf OUTPUT.tar.gz OUTPUT

        #  remove the viral DB
        rm -r virusbreakenddb_20210401
        rm -r BWA_index

    >>>

    output {
        File output1="~{sample_id}.virusbreakend.vcf.summary.tsv"
        File output2="~{sample_id}.virusbreakend.vcf"

    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(Human_Reference, "GB") + size(Virus_Reference, "GB") + size(fastq1, "GB")*6 + (disk)) + " HDD"
        docker: docker
        cpu: cpus
        memory: "100GB"
    }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

workflow VIRUSBreakend {
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
        # Reference Directories 
        #~~~~~~~~~~~~
        File Human_Reference
        File Virus_Reference


        #~~~~~~~~~~~~
        # general runtime settings
        #~~~~~~~~~~~~
        Int preemptible = 2
        String docker = "gridss/gridss:latest"
        Int disk = 100

        

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
    call RunVIRUSBreakend{
        input:
            fastq1 = left,
            fastq2 = right,


            Human_Reference = Human_Reference,
            Virus_Reference = Virus_Reference,
            
            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id,
            disk            = disk
    }
}
