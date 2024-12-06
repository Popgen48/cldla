///workflow to check input for at least three different input files///
///vcfsheet, mapsheet, outgroupsheet

workflow CHECK_INPUT {
    take:
        csvsheet

    main:
        Channel
            .fromPath(csvsheet)
            .splitCsv(header:true,sep:',')
            .map { create_individual_channel(it) }
            .set { meta_vcf_idx }

    emit:
        meta_vcf_idx
        }


def create_individual_channel(LinkedHashMap row){
        
        def meta = [:]
        def meta_vcf = []

        if (row.size() == 3 ){
            meta.id = row.chrom
            meta_vcf = [ meta, file(row.vcf), file(row.vcf_idx)  ]
        }
        else{
                meta.id = row.chrom
                meta_vcf = [ meta,file(row.file_path) ]
            }
        return meta_vcf
    }
