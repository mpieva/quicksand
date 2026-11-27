include { CHECK_DATABASE } from '../modules/local/check_database'

cyan = "\033[0;36m"
red = "\033[0;31m"
white = "\033[0m"
yellow = "\033[0;33m"

workflow setup {
    take: setup
    main:
        // write the commandline down
        log.info """
[quicksand]: Execution started: ${workflow.start.format('dd.MM.yyyy HH:mm')} ${cyan}

  ============================================================
  ==========================  =============================  =
  ==    ==  =  ==  ===   ===  =  ===   ====   ===  ========  =
  =  =  ==  =  ======  =  ==    ===  =  ==  =  ==     ===    =
  =  =  ==  =  ==  ==  =====   =====  =======  ==  =  ==  =  =
  ==    ==  =  ==  ==  =====    =====  ====    ==  =  ==  =  =
  ====  ==  =  ==  ==  =  ==  =  ==  =  ==  =  ==  =  ==  =  =
  ====  ===    ==  ===   ===  =  ===   ====    ==  =  ===    =
  ============================================================
  ${white}${workflow.manifest.description} ${cyan}~ Version ${workflow.manifest.version} ${white}

 --------------------------------------------------------------
        """

        //
        //
        // Additional Setup, Copy Files etc.
        //
        //
        def basedir = "quicksand_${workflow.manifest.version}"
        def start = workflow.start.format('yyyyMMdd_HHmmss')

        def commandFile = new File("${basedir}/nextflow/${start}_commands.txt")
        // create the nextflow dir if it doesnt exist yet
        if(!commandFile.getParentFile().exists()) {
            commandFile.getParentFile().mkdirs();
        }
        // write the pipeline run info to file
        commandFile << "# User\t${workflow.userName}\n"
        commandFile << "# Date/Time\t${workflow.start.format('yyyy-MM-dd HH:mm:ss')}\n"
        commandFile << "# Git Revision\t${workflow.commitId ?: 'local run'} (${workflow.revision ?: 'No version'})\n"
        commandFile << "${workflow.commandLine}\n"

        if( workflow.configFiles.size > 1 ){
            // configFiles[0] is the nextflow.config in the repository...
            // configFiles[1] is the one handed over with the -c flag
            def newconf = new File(workflow.configFiles[1] as String)
            def copyconf = new File("${basedir}/nextflow/${start}.config")

            copyconf << newconf.text
        }

        //
        // Check the database: do genomes and bedfiles match the taxonomy?
        // Mismatching genomes would be silently ignored in the mapping step
        //

        if( !params.rerun && params.db && params.genomes && params.bedfiles ){
            CHECK_DATABASE(
                file(params.db,       type:'dir', checkIfExists:true),
                file(params.genomes,  type:'dir', checkIfExists:true),
                file(params.bedfiles, type:'dir', checkIfExists:true),
            )

            CHECK_DATABASE.out.summary.subscribe{ summary ->
                def counts = summary.readLines()
                    .findAll{ it.trim() }
                    .collectEntries{ line ->
                        def fields = line.trim().split(/\s+/)
                        fields.size() > 1 && fields[1].isInteger() ? [(fields[0]): fields[1] as int] : [:]
                    }
                def issues = ['MISSING_GENOME', 'MISSING_BED'].findAll{ counts[it] }
                if( issues ){
                    log.warn "[quicksand]: ${yellow}Database check: ${issues.collect{ "${it}=${counts[it]}" }.join(', ')}. See stats/00_database_report.tsv${white}"
                } else {
                    log.info "[quicksand]: Database check: all genomes and bedfiles match the taxonomy"
                }
            }
        }
}