version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow run_fix_column {
    input {
        File vcf
        String sample_id 
    }
    
    call fix_column {
        input:
        vcf_file=vcf,
        sample_id=sample_id
    }
    output {
        File vcf_fixedcol=fix_column.vcf_fixedcol
    }
}

task fix_column {
    input {
        File vcf_file
        String sample_id

        Int memory = 2
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 10
        String docker = "python"
    }

    command {
    pip install bgzip

    echo """
import re
import sys
import bgzip
import gzip
print(sys.argv)
with gzip.open(sys.argv[1],'r+') as f:
    with open(sys.argv[2]+'_fixedcolumn.vcf.gz','wb') as raw2:
        with bgzip.BGZipWriter(raw2) as fout:
            for i, line in enumerate(f):
                original_string = line.decode('utf-8')
                if original_string[0] == '#':
                    fout.write(original_string.encode())
                else:
                    new_string = re.sub(
                        r'AS_FilterStatus=(.*?);', 
                        lambda x:'AS_FilterStatus=' + x.group(1).replace('|', '~').replace(',', '|').replace('~', ',') + ';', 
                        original_string)
                    to_print = new_string.split('\t')
                    to_print[7] = to_print[7].replace(' ','')
                    
                    #write the fixed string tab-separated to output file 
                    for k in range(len(to_print)-1):
                        fout.write((to_print[k] + '\t').encode())
                    fout.write(to_print[-1].encode())
    """ > script.py

    python script.py ${vcf_file} ${sample_id}
    }

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