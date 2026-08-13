
include { samplesheetToList } from 'plugin/nf-schema'

params {
    genome_version      : String  = "GRCh38"
    vcf_samplefile      : Path  = ""
    impute              : Boolean = false
    mocha_resources     : String = "/home/clarkb/images/mocha/resources"
    reference_fasta     : Path = "/home/clarkb/images/mocha/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    reference_fasta_fai : Path = "/home/clarkb/images/mocha/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
    beagle_resources    : Path = "${projectDir}/resources/GRCh38/beagle/"
}


process filterVCF {

    input:
    path vcf
    path reference
    path reference_index
    path segdups
    val chr

    output:
    tuple val("${chr}"), path("${vcf.simpleName}.filtered.vcf.gz" ), path("${vcf.simpleName}.filtered.vcf.gz.csi")


    script:
    """
    bcftools +mochatools --no-version -Ob ${vcf} -- -t GC -f ${reference} > vcf_withGC.bcf

    bcftools view --no-version -h vcf_withGC.bcf | sed 's/^\\(##FORMAT=<ID=AD,Number=\\)\\./\\1R/' | 
    bcftools reheader -h /dev/stdin vcf_withGC.bcf | 
    bcftools view -f 'PASS,.' |
    bcftools filter --no-version -Ou -e "FMT/DP<10 | FMT/GQ<20" --set-GT . | 
    bcftools view -T ^${segdups} -Ou |
    bcftools annotate --no-version -Ou -x ID,QUAL,^INFO/GC,^FMT/GT,^FMT/AD | 
    bcftools norm --no-version -Ou -m -any --keep-sum AD | 
    bcftools norm --no-version -o "${vcf.simpleName}.filtered.vcf.gz" -f "${reference}" --write-index -Oz -
    """
}

process phase {

    input:
    path vcf
    path index
    path reference_bref3
    path map_file
    val impute

    output:
    path "${vcf.simpleName}_BEAGLE_annotated.vcf.gz"

    script:
    """
    export JAVA_OPTS="-Xmx${task.memory.giga}g"


    java \$JAVA_OPTS -jar "${projectDir}/bin/beagle.27Feb25.75f.jar" gt=${vcf} \
        ref=${reference_bref3} \
        map=${map_file} \
        impute=${impute} \
        nthreads=${task.cpus} \
        out="${vcf.simpleName}_BEAGLE"

    tabix -p vcf "${vcf.simpleName}_BEAGLE.vcf.gz"

    bcftools annotate \
        --no-version \
        -o "${vcf.simpleName}_BEAGLE_annotated.vcf.gz" \
        -Oz --annotations "${vcf.simpleName}_BEAGLE.vcf.gz" \
        --columns -FMT/GT --write-index ${vcf}
    """

}
process sort {
    input: 
    path vcf
    path index

    output:
    path "${vcf.baseName}_sorted.bcf"

    script:
    """
    bcftools sort ${vcf} -b -o ${vcf.baseName}_sorted.bcf
    """


}

process mergeChromosomes {
    input:
    path vcf_list

    output:
    path "merged_phased.bcf", emit: "bcf"
    path  "merged_phased.bcf.csi", emit: "index"

    script:
    """
    for f in *.vcf.gz; do tabix -p vcf -f "\$f"; done

    bcftools concat --no-version -o merged_phased.bcf -Ob --write-index ${vcf_list}
    """
}



process mocha {

    input:
    path vcf
    path index
    path mocha_resources
    val genome_version

    output:
    path "mosaic_out.bcf", emit:"bcf"
    path "mosaic_stats.tsv", emit: "stats"
    path "mosaic_calls.tsv", emit: "calls"

    script:
    """

    bcftools +mocha \
            --genome ${genome_version} \
            --no-version \
            --output mosaic_out.bcf \
            --output-type b \
            --ucsc-bed mosaic_cnvs.bed \
            --calls mosaic_calls.tsv \
            --stats mosaic_stats.tsv \
            --write-index \
            --cnp ${mocha_resources}/cnps.bed \
            --mhc "chr6:27518932-33480487" \
            --kir "chr19:54071493-54992731" \
            --thr ${task.cpus} \
            ${vcf}
    """
}

process filterOutput {
    conda 'conda-forge::polars'

    input:
    path calls

    output:
    path "MCA_flagged.tsv"

    script:
    """
    filter_mocha.py ${calls}
    """

}
workflow {

main:
    bref3_ch = channel.fromPath("${params.beagle_resources}/bref3/*.bref3" )
        .map { file ->
                def chrom = (file.name =~ /chr(\w+)/)[0][1]  
                tuple(chrom, file)
        }

    map_ch = channel.fromPath("${params.beagle_resources}/chr_in_chrom_field/*.map" )
        .map { file ->
                def chrom = (file.name =~ /chrchr(\w+)/)[0][1]  
                tuple(chrom, file)
        }


    vcf_ch = channel.fromList(samplesheetToList(params.vcf_samplefile, "${projectDir}/assets/sampleFile_schema.json"))

    filterVCF(vcf_ch.map{it -> it[0]}, 
              params.reference_fasta, 
              params.reference_fasta_fai,
              file(params.mocha_resources).resolve("segdups.bed.gz"),
              vcf_ch.map{it -> it[1]})

    filtered_vcf_ch = filterVCF.out.join(bref3_ch).join(map_ch)
    
   

    phase_ch = phase(filtered_vcf_ch.map{it -> it[1]}, //vcf file
                     filtered_vcf_ch.map{it -> it[2]}, //index
                     filtered_vcf_ch.map{it -> it[3]}, //bref3 reference
                     filtered_vcf_ch.map{it -> it[4]}, //plink map file
                    params.impute)

    mergeChromosomes(phase_ch.collect())

    mocha(mergeChromosomes.out.bcf, 
          mergeChromosomes.out.index,  
          params.mocha_resources, 
          params.genome_version)

    filterOutput(mocha.out.calls)



publish:
    calls = filterOutput.out
    stats = mocha.out.stats
    raw_calls = mocha.out.calls
    
}

output {
    calls {

    }
    stats {

    }
    raw_calls{

    }

}
