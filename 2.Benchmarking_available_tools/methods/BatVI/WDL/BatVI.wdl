version 1.0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run BatVI
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
task RunBatVI {
    input {
        File fastq1
        File? fastq2

        # Human 
        File Human_BLAST_index
        File Human_BatIndex
        File Human_BWA_Index
        String Human_fasta

        # Virus 
        File Virus_BLAST_index
        File Virus_BatIndex
        File Virus_BWA_Index
        String Virus_fasta

        Int  insertion_length

        Int cpus
        Int preemptible
        String docker
        String sample_id
        Int disk
        String memory
    }

    command <<<
        set -ex


        #~~~~~~~~~~~~~~~~~~~~~~~~
        # Untar the references  
        #~~~~~~~~~~~~~~~~~~~~~~~~
        tar -xvf ~{Human_BLAST_index}
        tar -xvf ~{Human_BatIndex}
        tar -xvf ~{Human_BWA_Index}

        tar -xvf ~{Virus_BLAST_index}
        tar -xvf ~{Virus_BatIndex}
        tar -xvf ~{Virus_BWA_Index}


        #~~~~~~~~~~~~~~~~~~~~~~~~
        # Run BatVI
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

            /usr/local/src/adjust_read_names.py  --fastq1 ~{fastq1} --fastq2 ~{fastq2} --output .
      
            # Write configuration and filelist files
            /usr/local/src/BatVI/write_configuration_file.py \
                --fastq1 stripped_left.fq \
                --fastq2 stripped_right.fq \
                --Human_BLAST_Index ~{Human_BLAST_index} \
                --Human_BatIndex ~{Human_BatIndex} \
                --Human_BWA_Index ~{Human_BWA_Index} \
                --Human_fasta ~{Human_fasta} \
                --Virus_BLAST_Index ~{Virus_BLAST_index} \
                --Virus_BatIndex ~{Virus_BatIndex} \
                --Virus_BWA_Index ~{Virus_BWA_Index} \
                --Virus_fasta ~{Virus_fasta} \
                --insertion_length ~{insertion_length}
            
            #~~~~~~~~~~~~~~~~~~~~~~~
            # Run BatVI
            #~~~~~~~~~~~~~~~~~~~~~~~
            /usr/local/src/batvi1.03/call_integrations.sh `pwd` \
                --threads ~{cpus}  \
                2>&1 | tee output_log_subset.txt
        
        else 

            /usr/local/src/adjust_read_names.py  --fastq1 ~{fastq1} --fastq2 ~{fastq2} --output .
            
      
            # Write configuration and filelist files
            /usr/local/src/BatVI/write_configuration_file.py \
                --fastq1 stripped_left.fq \
                --fastq2 stripped_right.fq \
                --Human_BLAST_Index ~{Human_BLAST_index} \
                --Human_BatIndex ~{Human_BatIndex} \
                --Human_BWA_Index ~{Human_BWA_Index} \
                --Human_fasta ~{Human_fasta} \
                --Virus_BLAST_Index ~{Virus_BLAST_index} \
                --Virus_BatIndex ~{Virus_BatIndex} \
                --Virus_BWA_Index ~{Virus_BWA_Index} \
                --Virus_fasta ~{Virus_fasta} \
                --insertion_length ~{insertion_length}


            #~~~~~~~~~~~~~~~~~~~~~~~
            # Run BatVI
            #~~~~~~~~~~~~~~~~~~~~~~~
            /usr/local/src/batvi1.03/call_integrations.sh `pwd` \
                --threads ~{cpus}  \
                2>&1 | tee output_log_subset.txt
        fi



        mv batviconfig.txt ~{sample_id}.batviconfig.txt
        mv filelist.txt ~{sample_id}.filelist.txt
        mv final_hits.txt ~{sample_id}.final_hits.txt
        mv t.opt.subopt.cluster ~{sample_id}.t.opt.subopt.cluster
        mv clusterlist.opt.subopt.txt ~{sample_id}.clusterlist.opt.subopt.txt
        mv predictions.opt.subopt.txt ~{sample_id}.predictions.opt.subopt.txt



    >>>

    output {
        File batviconfig_txt="~{sample_id}.batviconfig.txt"
        File filelist_txt="~{sample_id}.filelist.txt"
        File final_hits_txt="~{sample_id}.final_hits.txt"
        File t_opt_subopt_cluster="~{sample_id}.t.opt.subopt.cluster"
        File clusterlist_opt_subopt_txt="~{sample_id}.clusterlist.opt.subopt.txt"
        File predictions_opt_subopt_txt="~{sample_id}.predictions.opt.subopt.txt"


    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(Human_BLAST_index, "GB") + size(Human_BatIndex, "GB") + size(Human_BWA_Index, "GB") + size(fastq1, "GB")*6 + disk) + " HDD"
        docker: docker
        cpu: cpus
        memory: memory + "GB"
    }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

workflow BatVI {
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
        # Human 
        File Human_BLAST_index
        File Human_BatIndex
        File Human_BWA_Index
        String Human_fasta

        # Virus 
        File Virus_BLAST_index
        File Virus_BatIndex
        File Virus_BWA_Index
        String Virus_fasta

        Int insertion_length 

        #~~~~~~~~~~~~
        # general runtime settings
        #~~~~~~~~~~~~
        Int preemptible = 2
        String docker = "trinityctat/batvi:devel"
        Int disk = 100
        String memory = "100"

        

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
    call RunBatVI{
        input:
            fastq1 = left,
            fastq2 = right,

            Human_BLAST_index = Human_BLAST_index,
            Human_BatIndex    = Human_BatIndex,
            Human_BWA_Index   = Human_BWA_Index,
            Human_fasta       = Human_fasta,

            Virus_BLAST_index = Virus_BLAST_index,
            Virus_BatIndex    = Virus_BatIndex,
            Virus_BWA_Index   = Virus_BWA_Index,
            Virus_fasta       = Virus_fasta,

            insertion_length  = insertion_length,
            
            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id,
            disk            = disk,
            memory          = memory
    }
}
