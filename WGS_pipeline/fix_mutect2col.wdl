# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow run_fix_column {
  call fix_column
}

task fix_column {
  File vcf
  String sample_id

  Int memory = 2
  Int boot_disk_size = 10
  Int num_threads = 1
  Int num_preempt = 5
  Int disk_space = 10
  String docker = "python"

  command <<<
    python <<CODE
    import re
    import gzip

    with gzip.open("${vcf}","r+") as f:
      with gzip.open("${sample_id}_fixedcolumn.vcf.gz","wb") as fout:
        for i, line in enumerate(f):
          original_string = line.decode('utf-8')
          if original_string[0] == "#":
            fout.write(original_string.encode())
          else:
            new_string = re.sub(
              r'AS_FilterStatus=(.*?);', 
              lambda x:"AS_FilterStatus=" + x.group(1).replace("|", "~").replace(",", "|").replace("~", ",") + ";", 
              original_string)
            fout.write(new_string.encode())
    CODE
  >>>

  output {
    File vcf_fixedcol="${sample_id}_fixedcolumn.vcf.gz"
  }

  runtime {
    docker: docker
    bootDiskSizeGb: "${boot_disk_size}"
    memory: "${memory}GB"
    disks: "local-disk ${disk_space} HDD"
    cpu: "${num_threads}"
    preemptible: "${num_preempt}"
  }

  meta {
    author: "Jeremie Kalfon"
  }
}