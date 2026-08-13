//TODO add to README how to add to resources folder
//      add documentation mocha resources

process getReferenceGenome {

    input:
    val chr

    output:
    tuple path("1kGP_high_coverage_Illumina.chr*.vcf.gz"), path("1kGP_high_coverage_Illumina.chr*.vcf.gz.tbi")


    script:
    """
    if [[ ${chr} == "X" ]]; then
        wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
        wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi"
    else
        wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"
    fi
    """

}

process fixPloidy {
    label "genomics"
    input:
    path chrX

    output:
    path "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel_diploid_only.v2.vcf.gz"

    script:
    """
    # Fix the chromosome X ploidy to phased diploid
    # Requires a ploidy.txt file containing 
    # space-separated CHROM,FROM,TO,SEX,PLOIDY 
    echo "chrX 1 156040895 M 2" > ploidy.txt
    bcftools +fixploidy \
        ${chrX}  -Ov -- -p ploidy.txt | 
        sed 's#0/0#0\\|0#g;s#1/1#1\\|1#g' | 
    bcftools view -Oz -o 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel_diploid_only.v2.vcf.gz
    """
}

process buildBref3{
    label "genomics"
    input:
    path vcf
    path index 

    output:
    path "*bref3"

    script:
    """

    java -jar ${projectDir}/bin/bref3.27Feb25.75f.jar ${vcf} > "${vcf.baseName}.bref3"

    """
}

process getMochaResources {

    output:
    path "GCA_000001405.15_GRCh38_no_alt_analysis_set.*"


    script:
    """
    wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | 
        gzip -d > GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    """
}

process index {
    input:
    path fasta

    output:
    path "*fai"

    script:
    """
    samtools faidx ${fasta}
    """


}
process getPlinkMap {

    output:
     path "chr_in_chrom_field/plink.chrchr*"

    script:
    """
    wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
    unzip plink.GRCh38.map.zip

    """
}


workflow {
    main:

    chr_all_ch = channel.of(1..22, 'X').flatten()
  

   // Build Bref3 reference panels
    ref_ch = getReferenceGenome(chr_all_ch)

    chr_ch = ref_ch.branch { vcf ->
        chrX : vcf =~ /^chrX/
        autosome : true

    }

    fixPloidy(chr_ch.chrX)

    all_chr_ch = chr_ch.autosome.mix(fixPloidy.out)


    buildBref3(all_chr_ch.map{it -> it[0]}, 
               all_chr_ch.map{it -> it[1]})
    //


    //Download Fasta
    getMochaResources()

    index(getMochaResources.out)

    //Get Plink Map files
    getPlinkMap()

    publish:
        plinks = getPlinkMap.out
        bref3 = buildBref3.out
        fasta = getMochaResources.out
        indices = index.out

}

output {
    plinks{
        mode 'copy'
        path "beagle"
    }
    bref3 {
        mode 'copy'
        path "beagle/bref3"
    }
    fasta {
        mode 'copy'
        path "mocha"
    }
    indices {
        mode 'copy'
        path 'mocha'
    }
}