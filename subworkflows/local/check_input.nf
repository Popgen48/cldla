///workflow to check input for at least three different input files///
///vcfsheet, mapsheet, outgroupsheet

workflow CHECK_INPUT{

    take:
        csvsheet

    main:
        Channel
            .fromPath(csvsheet)
            .splitCsv(sep:",")
            .map{ chrom, vcf, idx -> if(!file(vcf).exists() || !file(idx).exists()){ exit 1, "ERROR: Please check input vcfsheet, either vcf file or its index does not exit \
                -> ${vcf}" }else{tuple([id:chrom], file(vcf), file(idx))} }
            .set{ meta_vcf_idx }
       

    emit:
        meta_vcf_idx = meta_vcf_idx
}
