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
    --outdir                                Path to output results into (default: "results")
    --metadata_tsv                          Path to metadata TSV containing a column "name" with\n
                                            the filename of each FASTA input (without extension)\n
    --metadata_col                          Column in --metadata_tsv to use to annotate nodes/FASTA\n
                                            in the network visualisation\n
    --cpus								    Number of CPUs to assign to executor\n

	"""
}

process get_proteins {
    publishDir "results/1_orf_prediction", pattern: "*.faa"
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
    publishDir "results/2a_all_orfs", pattern: "all_orfs.faa"
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
    publishDir 'results/2b_annotations', pattern: "*.tsv"
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
    publishDir 'results/2b_annotations', pattern: "*.tsv"
    conda "$baseDir/conda_envs/blast.yml"
    input:
        path all_ORFs
    output:
        path "blast_vfs.tsv"
    script:
        """
        wget http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
        gunzip VFDB_setA_pro.fas.gz
        makeblastdb -in VFDB_setA_pro.fas -dbtype prot
        blastp -query $all_ORFs -db VFDB_setA_pro.fas -num_threads ${task.cpus} -outfmt 6 -max_target_seqs 1 -evalue 1e-25 > blast_vfs.tsv
        """
}


process run_proteinortho {
    publishDir 'results/3_proteinortho', pattern: "*.tsv"
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

process generate_graph_no_metadata {
    publishDir 'results/4_network_viz', pattern: "network.{pkl,png}"
    conda "$baseDir/conda_envs/graph.yml"
    input:
        tuple path(proteinortho), path(rgi_amr), path(vf_out6)
    output:
        tuple path("network.png"), path("network.pkl")
    script:
        """
        generate_graph.py --proteinortho_output $proteinortho --rgi_output $rgi_amr --vf_blast_output $vf_out6
        """
}

process generate_graph_metadata {
    publishDir 'results/4_network_viz', pattern: "network.{pkl,png}"
    conda "$baseDir/conda_envs/graph.yml"
    input:
        tuple path(proteinortho), path(rgi_amr), path(vf_out6), path(metadata_tsv)
        val(metadata_col)
    output:
        tuple path("network.png"), path("network.pkl")
    script:
        """
        generate_graph.py --proteinortho_output $proteinortho --rgi_output $rgi_amr --vf_blast_output $vf_out6 --seq_metadata $metadata_tsv --metadata_col $metadata_col
        """
}

process plot_interactive {
    publishDir 'results/4_network_viz', pattern: "network.html"
    conda "$baseDir/conda_envs/graph.yml"
    input:
        tuple path(network_png), path(network_pickle)
    output:
        path "network.html"
    script:
        """
        plot_interactive_network.py --networkx_pickle ${network_pickle}
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
    if (params.metadata_tsv && params.metadata_col ) {
        metadata = Channel.fromPath( params.metadata_tsv )
        graph_pickle = generate_graph_metadata( orthos
                                                .combine(amr)
                                                .combine(vf)
                                                .combine(metadata),
                                                params.metadata_col)
    } else {
        graph_pickle = generate_graph_no_metadata( orthos
                                                    .combine(amr)
                                                    .combine(vf) )
    }

    // use pyvis/visjs to create interactive version
    plot_interactive( graph_pickle )

}
