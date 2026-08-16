
include { samplesheetToList } from 'plugin/nf-schema'

params {
    genome_version      : String
    vcf_samplefile      : Path
    impute              : Boolean
    mocha_resources     : Path
    beagle_resources    : Path
    dataset_name        : String
}


process filterVCF {
    tag "Filtering variants from ${vcf.simpleName}"
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
    tag "Phasing ${vcf.simpleName}"
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


process mergeChromosomes {
    input:
    path vcf_list

    output:
    path "merged_phased.bcf", emit: "bcf"
    path "merged_phased.bcf.csi", emit: "index"

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
    path cnps
    val genome_version

    output:
    path "mosaic_out.bcf",   emit:"bcf"
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
            --cnp ${cnps} \
            --mhc "chr6:27518932-33480487" \
            --kir "chr19:54071493-54992731" \
            --thr ${task.cpus} \
            ${vcf}
    """
}

process filterOutput {
    input:
    path calls

    output:
    path "mosaicDB.tsv"

    script:
    """
    filter_mocha.py ${calls}
    """

}
workflow {

main:

    //Build channels from chromosome-level data for BEAGLE phasing
    bref3_ch = channel.fromPath("${params.beagle_resources}/bref3/*.bref3" )
        .map { file ->
                def chrom = (file.name =~ /(chr\w+)/)[0][1]  
                tuple(chrom, file)
        }

    map_ch = channel.fromPath("${params.beagle_resources}/chr_in_chrom_field/*.map" )
        .map { file ->
                def chrom = (file.name =~ /chr(chr\w+)/)[0][1]  
                tuple(chrom, file)
        }

    //Read samplesheet into chromosome-level channel
    vcf_ch = channel.fromList(samplesheetToList(params.vcf_samplefile, "${projectDir}/assets/sampleFile_schema.json"))

    //Filter the VCF as recommended by MoChA
    filterVCF(vcf_ch.map{it -> it[0]},                                                                          //chromosome num
              file(params.mocha_resources).resolve("GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"),          //fasta reference for GC mapping
              file(params.mocha_resources).resolve("GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"),      //index
              file(params.mocha_resources).resolve("segdups.bed.gz"),                                           //segdup regions
              vcf_ch.map{it -> it[1]})                                                                          // the vcf file

    // Join BEAGLE chromosome specific resource into one tuple
    filtered_vcf_ch = filterVCF.out.join(bref3_ch).join(map_ch)
    
   
    // Perform BEAGLE phasing
    phase_ch = phase(filtered_vcf_ch.map{it -> it[1]}, //vcf file
                     filtered_vcf_ch.map{it -> it[2]}, //index
                     filtered_vcf_ch.map{it -> it[3]}, //bref3 reference
                     filtered_vcf_ch.map{it -> it[4]}, //plink map file
                    params.impute)

    //Pull all phased chromosome files into one list and merge them into one bcf
    mergeChromosomes(phase_ch.collect())

    //Run MoChA
    mocha(mergeChromosomes.out.bcf, 
          mergeChromosomes.out.index,  
          file(params.mocha_resources).resolve("cnps.bed"), 
          params.genome_version)


    //Make filtered output for high confidence somatic mutations with more readable column names
    filterOutput(mocha.out.calls)



publish:
    calls = filterOutput.out
    stats = mocha.out.stats
    raw_calls = mocha.out.calls
    
}

output {
    calls {
        path "${params.dataset_name}"
    }
    stats {
        path "${params.dataset_name}"
    }
    raw_calls{
        path "${params.dataset_name}"
    }

}
