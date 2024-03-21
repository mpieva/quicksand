workflow write_reports {
    take: ch_final
    take: ch_versions
    main:

    def basedir = "quicksand_${workflow.manifest.version}"

    //
    //
    // Write Reports
    //
    //

    // write the reports to file...
    ch_versions.unique().collectFile(name: 'pipeline_versions.yml', storeDir:"${basedir}/")

    // calculate proportion mapped, duplication rate, expected breadth
    ch_final
    .map{ meta -> meta+[
            "ProportionMapped": ( (meta.ReadsExtracted as int)==0 || (meta.ReadsMapped as int)== 0) ? 0 : ((meta.ReadsMapped as int)/(meta.ReadsExtracted as int)).trunc(4),
            "DuplicationRate":  ( (meta.ReadsDeduped as int)  ==0 || (meta.ReadsMapped as int)== 0) ? 0 : ((meta.ReadsMapped as int)/(meta.ReadsDeduped   as int)).trunc(4),
            "ReadsFinal": meta.ReadsBedfiltered == '-' ? meta.ReadsDeduped as int : meta.ReadsBedfiltered as int,
        ]
    }
    .unique{meta -> meta.RG+meta.Species+meta.Reference}
    .set{ ch_final }

    //
    // calculate the FamPercentage
    //
    //         best/fixed
    //             |
    //    maximum value (by family and rg)
    //    (no difference for best)
    //             |
    //    sum the family reads by RG = total_rg
    //             |
    //    add to ch_final by RG
    //             |
    //    calculate FamPercent by Species
    //

    ch_final
    .map{meta ->
        [
            [meta.RG, meta.Family],
            meta.ReadsFinal,
            meta
        ]
    }
    .groupTuple()
    .map{key, vals, metas ->
        [key, ['FamilyMeanReads':vals.max()]]
    } // now set key to RG
    .map{ key, val ->
        [key[0], val.FamilyMeanReads]
    } //and sum up
    .groupTuple()
    .map{ key, val ->
        [key, val.sum()]
    }
    .set{ total_rg }

    // and combine with the main channel
    ch_final
    .map{ meta  -> [ meta.RG, meta ] } //prepare id as key
    .combine(total_rg, by:0)
    .map{rg, meta, total_rg -> // and calculate
        meta+['FamPercentage': total_rg==0 ? 0 : (meta.ReadsFinal * 100 / total_rg).trunc(2)]
    }
    .set{ch_final}

    //
    // Write now all the data to files!
    //

    header_map = [
    'tax'    : 'Order\tFamily\tSpecies\tReference',
    'split'  : 'ReadsRaw\tReadsFiltered\tReadsLengthfiltered',
    'kraken' : 'SpeciesKmers\tKmerCoverage\tKmerDupRate',
    'deam'   : 'Ancientness\tReadsDeam(1term)\tReadsDeam(3term)\tDeam5(95ci)\tDeam3(95ci)\tDeam5Cond(95ci)\tDeam3Cond(95ci)',
    'extract': 'ExtractLVL\tReadsExtracted',
    'map'    : 'ReadsMapped\tProportionMapped',
    'dedup'  : 'ReadsDeduped\tDuplicationRate\tCoveredBP',
    'bed'    : 'ReadsBedfiltered\tPostBedCoveredBP',
    'frags'  : 'MeanFragmentLength\tMeanFragmentLength(3term)',
    'breadth': 'Coverage\tBreadth\tExpectedBreadth\tProportionExpectedBreadth'
    ]

    def getVals = {String header, meta, res=[] ->
        header.split('\t').each{res << meta[it]}
        res.join('\t')
    }

    ch_final
    .collectFile( name:"final_report.tsv",
        seed:[
        'RG',
        header_map['split'],
        header_map['kraken'],
        header_map['extract'],
        header_map['tax'],
        header_map['map'],
        header_map['dedup'],
        header_map['bed'],
        'FamPercentage',
        header_map['deam'],
        header_map['frags'],
        header_map['breadth']
        ].join('\t'), storeDir:"${basedir}/", newLine:true, sort:true
    ){[
        it.RG,
        getVals(header_map['split'],   it),
        getVals(header_map['kraken'],  it),
        getVals(header_map['extract'], it),
        getVals(header_map['tax'],     it),
        getVals(header_map['map'],     it),
        getVals(header_map['dedup'],   it),
        getVals(header_map['bed'],     it),
        it.FamPercentage,
        getVals(header_map['deam'],    it),
        getVals(header_map['frags'],   it),
        getVals(header_map['breadth'], it)
        ].join('\t')
    }
    .subscribe {
        println "[quicksand]: Summary reports saved"
    }

    ch_final
    .filter{ it.FamPercentage >= params.reportfilter_percentage  }
    .filter{ it.ProportionExpectedBreadth >= params.reportfilter_breadth }
    .collectFile( name:"filtered_report_${params.reportfilter_percentage}p_${params.reportfilter_breadth}b.tsv",
        seed:[
        'RG',
        header_map['split'],
        header_map['kraken'],
        header_map['extract'],
        header_map['tax'],
        header_map['map'],
        header_map['dedup'],
        header_map['bed'],
        'FamPercentage',
        header_map['deam'],
        header_map['frags'],
        header_map['breadth']
        ].join('\t'), storeDir:"${basedir}/", newLine:true, sort:true
    ){[
        it.RG,
        getVals(header_map['split'],   it),
        getVals(header_map['kraken'],  it),
        getVals(header_map['extract'], it),
        getVals(header_map['tax'],     it),
        getVals(header_map['map'],     it),
        getVals(header_map['dedup'],   it),
        getVals(header_map['bed'],     it),
        it.FamPercentage,
        getVals(header_map['deam'],    it),
        getVals(header_map['frags'],   it),
        getVals(header_map['breadth'], it)
        ].join('\t')
    }
    .subscribe {
        println "[quicksand]: Filtered report saved"
    }

    ch_final
    .collectFile(
        storeDir: "${basedir}/stats", newLine:true,
        seed:[
        header_map['tax'],
        header_map['deam']
        ].join('\t')
    ){[
        "${it.RG}_04_deamination.tsv", [
            getVals(header_map['tax'],  it),
            getVals(header_map['deam'], it)
        ].join('\t')
        ]
    }

    ch_final
    .collectFile(
        storeDir: "${basedir}/stats", newLine:true,
        seed: [
        header_map['tax'],
        header_map['bed']
        ].join('\t')
    ){[
        "${it.RG}_03_bedfiltered.tsv",
            [
            getVals(header_map['tax'], it),
            getVals(header_map['bed'], it)
            ].join('\t')
        ]
    }

    ch_final
    .collectFile(
        storeDir: "${basedir}/stats", newLine:true,
        seed: [
        header_map['tax'],
        header_map['dedup'],
        header_map['breadth']
        ].join('\t')
    ){[
        "${it.RG}_02_deduped.tsv",
            [
            getVals(header_map['tax'],   it),
            getVals(header_map['dedup'], it),
            getVals(header_map['breadth'], it)
            ].join('\t')
        ]
    }

    ch_final
    .collectFile(
        storeDir: "${basedir}/stats", newLine:true,
        seed: [
        header_map['tax'],
        header_map['map']
        ].join('\t')
    ){[
        "${it.RG}_01_mapped.tsv",
            [
            getVals(header_map['tax'], it),
            getVals(header_map['map'], it)
            ].join('\t')
        ]
    }

    ch_final
    .unique{ meta -> meta.RG+meta.Taxon }
    .collectFile(storeDir: "${basedir}/stats", seed:'Taxon\tReadsExtracted', newLine:true) {meta ->
        [ "${meta.RG}_00_extracted.tsv", "${meta.Taxon}\t${meta.ReadsExtracted}"]
    }

    ch_final
    .unique{it.RG}
    .collectFile(
        storeDir: "${basedir}/stats", newLine:true,
        seed: [
        'RG',
        header_map['split'],
        ].join('\t')
    ){[
        "splitcounts.tsv",
            [
            it.RG,
            getVals(header_map['split'], it),
            ].join('\t')
        ]
    }
}
