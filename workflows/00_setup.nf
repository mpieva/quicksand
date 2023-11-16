workflow setup {
    take: setup
    main:
        // write the commandline down

        def start = workflow.start.format('yyyyMMdd_HHmmss')

        def commandFile = new File("nextflow/${start}_commands.txt")
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
            def copyconf = new File("nextflow/${start}.config")

            copyconf << newconf.text
        }
}