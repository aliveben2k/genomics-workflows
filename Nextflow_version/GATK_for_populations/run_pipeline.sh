#!/usr/bin/env bash

set -euo pipefail

pipeline_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base_workdir="${PWD}/work"
base_outdir="${PWD}/results"
run_serial=""
run_name=""
resume_target=""
is_resume=false
serial_was_supplied=false
forward_args=()

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_value() {
    [[ $# -ge 2 && -n "$2" ]] || die "Option $1 requires a value."
}

print_help() {
    cat <<'EOF'
nf-core/gatkpopulation launcher

USAGE
  ./run_pipeline.sh [Nextflow options] [pipeline parameters]

The launcher appends a four-character run serial to the base work and output
directories, stores nextflow.log with the work directory, and forwards all
unrecognized arguments to Nextflow.

HELP AND CONFIGURATION
  -h, --help, --help-full
      Print this help and exit without creating a run.

  -profile PROFILE[,PROFILE...]
      Select one executor profile and one software profile. Add gpu to use
      Parabricks. Examples:
        -profile slurm,conda
        -profile sge,singularity
        -profile uge,conda
        -profile pbs,conda
        -profile pbspro,singularity
        -profile slurm,singularity,gpu

  -params-file FILE
      YAML or JSON analysis parameters. Start from:
        assets/analysis.template.yml

  -c FILE
      Nextflow process-resource configuration. Start from:
        assets/resources.template.config

LAUNCHER PATHS AND RUN IDENTITY
  -work-dir DIR
      Base Nextflow work directory. Default: ./work
      Serial K8Z2 produces DIR_K8Z2.

  --outdir DIR
      Base published-results directory. Default: ./results
      Serial K8Z2 produces DIR_K8Z2.

  --run_serial CODE, --run-serial CODE
      Optional four-character run serial matching [A-Z][A-Z0-9]{3}.
      An unused serial is generated when omitted.

  -name NAME
      Nextflow run name. The default is gatkpopulation_SERIAL.

  -resume NAME
      Resume the named run. The launcher recovers its serial from the marker
      stored in the serial-specific work directory.

INPUT AND STAGE CONTROL
  --input FILE
      Stage-specific CSV manifest.

  --start_stage N
      First stage to run. Default: 0.

  --stop_stage N
      Last stage to run. Default: 6.

  Stage  Operation                         Input CSV columns
  0      Reference preparation            sample,fastq_1,fastq_2
  1      Trimming and FastQC               sample,fastq_1,fastq_2
  2      BAM preprocessing                 sample,fastq_1,fastq_2
  3      Per-sample GVCF                   sample,bam,bai
  4      Per-contig GenomicsDB             sample,gvcf,tbi
  5      Per-contig joint genotyping       contig,genomicsdb
  6      VCF filtering                     contig,vcf,tbi

REFERENCE OPTIONS
  --fasta FILE
      Uncompressed reference FASTA.

  --fasta_fai FILE
      Existing FASTA index. Inferred as FASTA.fai when starting after stage 0.

  --fasta_dict FILE
      Existing GATK sequence dictionary. Inferred beside the FASTA.

  --bwa_index DIR
      Existing traditional BWA index directory. Required when starting at
      stage 1 or 2 and running stage 2.

  --known_sites FILE[,FILE...]
      Indexed known-sites VCFs for BQSR. BQSR is skipped when omitted.
      A YAML list in the parameter file is recommended for multiple files.

  --intervals FILE
      Optional intervals for BQSR and per-sample GVCF calling.

TRIMMING OPTIONS
  --trimmer fastp|trimmomatic
      Trimmer to run before FastQC. Default: fastp.

  --adapter_fasta FILE
      Optional adapter FASTA supplied to fastp.

  --trimmomatic_args STRING
      Ordered Trimmomatic steps.
      Default: SLIDINGWINDOW:4:20 MINLEN:36

SEQUENCING TYPE
  --sequencing_type wgs|wes|rad|amp
      Input design. Default: wgs.
      wgs/wes run duplicate marking and produce *_aln_sort_MD.bam.
      rad/amp retain duplicates and produce *_aln_sort.bam.
      rad covers RAD-seq, GBS, ddRAD, and related reduced-representation data.

  --wgs, --wes, --rad, --amp
      Convenience aliases for --sequencing_type. Use only one.

GPU OPTIONS
  Add the gpu profile with Docker, Singularity, or Apptainer to use Parabricks.

  --parabricks_sif FILE
      Local Parabricks .sif for Singularity or Apptainer.

  --gpu_count N
      GPUs per Parabricks task. Default: 4.

  --gpu_fallback true|false
      Rerun unsuccessful GPU tasks with CPU tools. Default: true.

  --gpu_max_submit_await DURATION
      Maximum GPU queue wait before fallback. Default: 30 min.

JOINT GENOTYPING
  --all_sites
      Include nonvariant loci. Per-contig raw names become
      CONTIG.all_sites.raw.vcf.gz.

  --subset_samples FILE
      One-column CSV with header "sample". SelectVariants extracts those
      samples from each GenomicsDB under 05.1-subset before GenotypeGVCFs.

VCF FILTERING
  --filter_vcf biallele|monobi
      Named stage-6 preset. Default: biallele.
      biallele retains quality-filtered biallelic SNPs.
      monobi retains monomorphic sites and quality-filtered biallelic SNPs,
      while excluding indels. Use monobi with --all_sites unless stage-6 input
      already contains monomorphic records.

  --bcftools_view_args STRING
      Advanced override for the selected preset. The pipeline always appends
      compressed output and tabix indexing options. An empty string retains
      every input record.

  --gather_vcfs 5|6
      5: gather raw contig VCFs into 05-raw_vcf/all.raw.vcf.gz, then filter
         only the gathered VCF if stage 6 runs.
      6: filter contigs first, then gather into
         06-filtered_vcf/all.filtered.vcf.gz.

SCHEDULER OPTIONS
  Slurm profile:
    --slurm_queue NAME
    --slurm_account NAME
    --slurm_qos NAME
    --slurm_queue_size N

  SGE profile:
    --sge_queue NAME
    --sge_project NAME
    --sge_penv NAME
    --sge_cluster_options STRING
    --sge_gpu_cluster_options STRING
    --sge_queue_size N

  Univa Grid Engine profile:
    --uge_queue NAME
    --uge_project NAME
    --uge_penv NAME
    --uge_cluster_options STRING
    --uge_gpu_cluster_options STRING
    --uge_queue_size N

  PBS/Torque profile:
    --pbs_queue NAME
    --pbs_account NAME
    --pbs_cluster_options STRING
    --pbs_gpu_cluster_options STRING
    --pbs_queue_size N

  PBS Pro profile:
    --pbspro_queue NAME
    --pbspro_account NAME
    --pbspro_cluster_options STRING
    --pbspro_gpu_cluster_options STRING
    --pbspro_queue_size N

  GPU cluster-option strings can use {gpu_count}, which is replaced with the
  value of --gpu_count.

REPORTING AND OUTPUT
  --publish_dir_mode copy|copyNoFollow|link|move|rellink|symlink
      Nextflow publishing mode. Default: copy.

  --multiqc_title STRING
  --multiqc_config FILE
  --multiqc_logo FILE
  --multiqc_methods_description FILE
      Optional MultiQC customization.

EXAMPLES
  CPU WGS on Slurm:
    ./run_pipeline.sh \
      -profile slurm,conda \
      -params-file /PATH/TO/analysis.yml \
      -c /PATH/TO/resources.config \
      -work-dir /PATH/TO/work/gatkpopulation \
      --outdir /PATH/TO/results/gatkpopulation

  GPU RAD-seq with Singularity:
    ./run_pipeline.sh \
      -profile slurm,singularity,gpu \
      --rad \
      --parabricks_sif /PATH/TO/parabricks.sif \
      --input /PATH/TO/samplesheet.csv \
      --fasta /PATH/TO/reference.fasta \
      -work-dir /PATH/TO/work/gatkpopulation \
      --outdir /PATH/TO/results/gatkpopulation

  Resume:
    ./run_pipeline.sh \
      -resume gatkpopulation_K8Z2 \
      -profile slurm,conda \
      -work-dir /PATH/TO/work/gatkpopulation \
      --outdir /PATH/TO/results/gatkpopulation

See README.md, docs/usage.md, nextflow_schema.json, and the templates under
assets/ for the complete configuration reference.
EOF
}

random_serial() {
    local letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    local alphanumeric='ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
    printf '%s%s%s%s' \
        "${letters:RANDOM%26:1}" \
        "${alphanumeric:RANDOM%36:1}" \
        "${alphanumeric:RANDOM%36:1}" \
        "${alphanumeric:RANDOM%36:1}"
}

marker_value() {
    local marker="$1"
    local key="$2"
    awk -F= -v key="$key" '$1 == key { sub(/^[^=]*=/, ""); print; exit }' "$marker"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help|--help-full|--help_full)
            print_help
            exit 0
            ;;
        --outdir)
            require_value "$@"
            base_outdir="$2"
            shift 2
            ;;
        --outdir=*)
            base_outdir="${1#*=}"
            shift
            ;;
        -work-dir)
            require_value "$@"
            base_workdir="$2"
            shift 2
            ;;
        -work-dir=*)
            base_workdir="${1#*=}"
            shift
            ;;
        --run_serial|--run-serial)
            require_value "$@"
            run_serial="$2"
            serial_was_supplied=true
            shift 2
            ;;
        --run_serial=*|--run-serial=*)
            run_serial="${1#*=}"
            serial_was_supplied=true
            shift
            ;;
        -name)
            require_value "$@"
            run_name="$2"
            forward_args+=("$1" "$2")
            shift 2
            ;;
        -resume)
            is_resume=true
            forward_args+=("$1")
            if [[ $# -gt 1 && "$2" != -* ]]; then
                resume_target="$2"
                forward_args+=("$2")
                shift 2
            else
                shift
            fi
            ;;
        *)
            forward_args+=("$1")
            shift
            ;;
    esac
done

if [[ "$is_resume" == true && -z "$resume_target" && -z "$run_name" ]]; then
    die "Use -resume <run_name> so the correct concurrent session is selected."
fi

if [[ -z "$run_serial" && "$is_resume" == true ]]; then
    lookup_name="${resume_target:-$run_name}"
    [[ -n "$lookup_name" ]] || die \
        "Resume without --run_serial requires an explicit run name: -resume <run_name>."

    shopt -s nullglob
    matching_markers=()
    for marker in "${base_workdir}"_*/.gatkpopulation-run; do
        if [[ "$(marker_value "$marker" run_name)" == "$lookup_name" ]]; then
            matching_markers+=("$marker")
        fi
    done
    shopt -u nullglob

    [[ ${#matching_markers[@]} -eq 1 ]] || die \
        "Expected one stored serial for run '$lookup_name' under ${base_workdir}_*, found ${#matching_markers[@]}."
    run_serial="$(marker_value "${matching_markers[0]}" run_serial)"
fi

if [[ -z "$run_serial" ]]; then
    for _attempt in $(seq 1 1000); do
        candidate="$(random_serial)"
        if [[ ! -e "${base_workdir}_${candidate}" && ! -e "${base_outdir}_${candidate}" ]]; then
            run_serial="$candidate"
            break
        fi
    done
    [[ -n "$run_serial" ]] || die "Unable to generate an unused run serial after 1000 attempts."
fi

[[ "$run_serial" =~ ^[A-Z][A-Z0-9]{3}$ ]] || die \
    "Run serial '$run_serial' must match ^[A-Z][A-Z0-9]{3}$."

scoped_workdir="${base_workdir}_${run_serial}"
scoped_outdir="${base_outdir}_${run_serial}"

if [[ "$is_resume" == false && "$serial_was_supplied" == true ]] &&
   [[ -e "$scoped_workdir" || -e "$scoped_outdir" ]]; then
    die "Serial '$run_serial' already has a work or output directory. Use -resume or choose another serial."
fi

if [[ "$is_resume" == false && -z "$run_name" ]]; then
    run_name="gatkpopulation_${run_serial}"
    forward_args+=("-name" "$run_name")
fi

mkdir -p "$scoped_workdir" "$scoped_outdir"

marker_run_name="${run_name:-$resume_target}"
if [[ -n "$marker_run_name" ]]; then
    {
        printf 'run_serial=%s\n' "$run_serial"
        printf 'run_name=%s\n' "$marker_run_name"
        printf 'base_workdir=%s\n' "$base_workdir"
        printf 'base_outdir=%s\n' "$base_outdir"
    } > "${scoped_workdir}/.gatkpopulation-run"
fi

printf 'Run serial: %s\n' "$run_serial"
printf 'Work directory: %s\n' "$scoped_workdir"
printf 'Output directory: %s\n' "$scoped_outdir"

exec nextflow \
    -log "${scoped_workdir}/nextflow.log" \
    run "$pipeline_dir" \
    "${forward_args[@]}" \
    -work-dir "$scoped_workdir" \
    --outdir "$scoped_outdir" \
    --run_serial "$run_serial"
