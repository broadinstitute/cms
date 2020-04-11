version 1.0

import "tasks_simulation.wdl" as sims

workflow run_sims {
    input {
      Array[File] paramFiles
      File recombFile
      Int nreps
      String       cms_docker = "quay.io/ilya_broad/cms"
    }

    scatter(paramFile in paramFiles) {
        scatter(rep in range(nreps)) {
            call sims.cosi2_run_one_sim {
                input:
                   paramFile = paramFile, recombFile=recombFile, simId=basename(paramFile, ".par")+"_"+rep, cms_docker=cms_docker
            }
        }
    }
}


