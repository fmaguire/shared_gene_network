#!/usr/bin/env nextflow 
nextflow.enable.dsl=2

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
        blastp -query $all_ORFs -db VFDB_setB_pro.fas -num_threads ${task.cpus} -outfmt 6 -max_target_seqs 1 -evalue 1e-10 > blast_vfs.tsv
        """
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

    input_ch = Channel.fromPath( params.input_sequences + "/*.{fa,fasta,fna}" )
    get_proteins( input_ch )
    proteins = get_proteins.out.collect()
    
    combine_all_ORFs( proteins )      
    ORFs = combine_all_ORFs.out

    amr = annotate_amr( ORFs )
    vf =  annotate_vf( ORFs )

    orthos = run_proteinortho( proteins )
    
    generate_graph( orthos.combine(amr).combine(vf) )


}
