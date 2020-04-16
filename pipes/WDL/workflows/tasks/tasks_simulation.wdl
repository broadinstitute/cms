version 1.0

task cosi2_run_one_sim {
  input {
    File         paramFile
    File         recombFile
    String       simId
    Int          maxAttempts = 10000000
    String       cms_docker = "quay.io/ilya_broad/cms"
  }

  command {
    
    grep -v "recomb_file" "${paramFile}" > ${paramFile}.fixed.par
    echo "recomb_file ${recombFile}" >> ${paramFile}.fixed.par
    env COSI_NEWSIM=1 COSI_MAXATTEMPTS=${maxAttempts} coalescent -p ${paramFile}.fixed.par --genmapRandomRegions --drop-singletons .25 --output-gen-map --tped "${simId}.tped"
  }

  output {
    Array[File]        tped = glob("*.tped")
    String      cms_version = "cms_version_unknown"
  }
  runtime {
    docker: "${cms_docker}"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}


