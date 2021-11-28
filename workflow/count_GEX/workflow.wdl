version 1.0

struct RuntimeAttr {
  Int? cpu
  Int? memory_gb
  String? docker
  String? queue
}

workflow cellranger_count_GEX {
  input {
    Array[File] fq1s
    Array[File] fq2s
    String sample_name
    Int expect_cells = 3000
    String transcriptome

    RuntimeAttr runtime_attr_count
  }
  call count {
    input:
      fq1s = fq1s,
      fq2s = fq2s,
      sample_name = sample_name,
      expect_cells = expect_cells,
      transcriptome = transcriptome,

      runtime_attr_override = runtime_attr_count
  }
  output {
    File report = count.out_report
    File bam = count.out_bam
    File h5ad = count.out_h5ad
  }
}

task count {
  input {
    Array[File] fq1s
    Array[File] fq2s
    String sample_name
    Int expect_cells
    String transcriptome

    RuntimeAttr runtime_attr_override
  }
  command {
    set -euo pipefail
    python3 /opt/10x/raw_to_tenx.py --fq1s ~{sep=" --fq1s " fq1s} --fq2s ~{sep=" --fq2s " fq2s} --sample_name "~{sample_name}"
    cellranger count --id=~{sample_name} --fastqs=rawdata --sample=~{sample_name} --expect-cells=~{expect_cells} --transcriptome=~{transcriptome} --localmem=~{memory_gb} --nopreflight --disable-ui
    python3 /opt/10x/tenx_to_h5ad.py --h5 "~{sample_name}/outs/filtered_feature_bc_matrix.h5"
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 32
  }
  Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
  Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
  output {
    File out_report = "~{sample_name}/outs/web_summary.html"
    File out_h5ad = "~{sample_name}/outs/filtered_feature_bc_matrix.h5ad"
    File out_bam = "~{sample_name}/outs/possorted_genome_bam.bam"
  }
}
