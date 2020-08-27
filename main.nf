#!/usr/bin/env nextflow 
nextflow.enable.dsl=2

def helpMessage() {
    log.info """

\nUsage:
The typical command for running the pipeline is as follows:
nextflow run fmaguire/shared_gene_network --input_sequences data
Mandatory arguments:
    --input_sequences                       Path to input data: folder containing FASTA formatted genomes/contigs

Optional arguments:
    --cpus								    Number of CPUs to assign to executor\n

	"""
}

process get_proteins {
    publishDir "results/1_orf_prediction", pattern: "*.faa", mode: 'copy'
    conda "$baseDir/conda_envs/prodigal.yml"
    input:
        path fasta 
    output:
        path "${fasta.simpleName}.faa" 
    script:
        """
        prodigal -p meta -a ${fasta.simpleName}.faa -i ${fasta}
        """
}

process combine_all_ORFs {
    publishDir "results/2a_all_orfs", pattern: "all_orfs.faa", mode: 'copy'
    input:
        path all_fasta_aa
    output:
        path "all_orfs.faa"
    script:
        """
        cat $all_fasta_aa > all_orfs.faa
        """
}


process annotate_amr {
    publishDir 'results/2b_annotations', pattern: "*.tsv", mode: 'copy'
    conda "$baseDir/conda_envs/rgi.yml"
    input:
        path all_ORFs
    output:
        path "rgi_amr.tsv"
    script:
        """
        tr -d '*' < $all_ORFs > all_orfs_clean.faa
        rgi main -i all_orfs_clean.faa -o rgi_amr -t protein -n ${task.cpus} 
        mv rgi_amr.txt rgi_amr.tsv
        """
}


process annotate_vf {
    publishDir 'results/2b_annotations', pattern: "*.tsv", mode: 'copy'
    conda "$baseDir/conda_envs/blast.yml"
    input:
        path all_ORFs
    output:
        path "blast_vfs.tsv"
    script:
        """
        wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz 
        gunzip VFDB_setB_pro.fas.gz
        makeblastdb -in VFDB_setB_pro.fas -dbtype prot
        blastp -query $all_ORFs -db VFDB_setB_pro.fas -num_threads ${task.cpus} -outfmt 6 -max_target_seqs 1 -evalue 1e-25 > blast_vfs.tsv
        """
        //sed -i 's|) (|)_(|' VFDB_setB_pro.fas
}


process run_proteinortho {
    publishDir 'results/3_proteinortho', pattern: "*.tsv", mode: 'copy'
    conda "$baseDir/conda_envs/proteinortho.yml"
    input:
        path all_aa
    output: 
        path "network.proteinortho.tsv"
    script:
        """
        proteinortho -project=network -singles -cpus=${task.cpus} -silent -clean $all_aa
        """
}


process generate_graph {
    publishDir 'results/4_network_viz', pattern: "network.png", mode: 'copy'
    conda "$baseDir/conda_envs/graph.yml"
    input:
        tuple path(proteinortho), path(rgi_amr), path(vf_out6)
    output:
        path "network.png"
    script:
        """
        generate_graph.py --proteinortho_output $proteinortho --rgi_output $rgi_amr --vf_blast_output $vf_out6
        """
}


workflow {

	if (params.help) {
		helpMessage()
		exit 0
	}
    
    if (!params.input_sequences) {
        println("\nMust provide folder of FASTA contigs/genomes/plasmids as --input_sequences\n")
        helpMessage()
        exit 1
    }
    
    // Detect ORFs and translate proteins using prodigal
    input_ch = Channel.fromPath( params.input_sequences + "/*.{fa,fasta,fna}" )
                .ifEmpty { error "\nCannot find any FASTA formatted sequences in ${params.input_sequences}\n" }
    
    get_proteins( input_ch )
    proteins = get_proteins.out.collect()
    
    // Collate all ORFs and annotate them for AMR (using CARD+RGI)
    // and VFs (using VFDB+BLASTP)
    ORFs = combine_all_ORFs( proteins )      
    amr = annotate_amr( ORFs )
    vf =  annotate_vf( ORFs )
    
    // Predict orthologues 
    orthos = run_proteinortho( proteins )
    
    // Use networkx script to generate shared gene network and highlight with
    // amr and vf annotations
    generate_graph( orthos.combine(amr).combine(vf) )

}
