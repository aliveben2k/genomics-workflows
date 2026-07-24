#!/usr/bin/env python3

import sys
import argparse
import io
import re
import subprocess
import tempfile
from contextlib import ExitStack
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pysam

def read_vcf_header(vcf_file):
    header_lines = []
    ids = None
    contig_lengths = {}
    try:
        with subprocess.Popen(
            ["bcftools", "view", "-h", vcf_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ) as proc:
            assert proc.stdout is not None
            for line in proc.stdout:
                if line.startswith("##contig="):
                    m = re.search(r"ID=([^,\s>]+)", line)
                    if m:
                        contig_id = m.group(1).strip()
                        lm = re.search(r"length\s*=\s*(\d+)", line, re.IGNORECASE)
                        if lm:
                            contig_lengths[contig_id] = int(lm.group(1))
                if line.startswith("#CHROM"):
                    ids = line.rstrip("\r\n").split("\t")
                    break
                if line.startswith("#"):
                    header_lines.append(line.rstrip("\r\n"))
            return_code = proc.wait()
            if return_code != 0:
                stderr = proc.stderr.read().strip() if proc.stderr else ""
                message = f"bcftools failed with exit code {return_code}"
                if stderr:
                    message += f": {stderr}"
                raise RuntimeError(message)
    except OSError as e:
        raise RuntimeError(f"Cannot read VCF header from {vcf_file}: {e}") from e
    if not ids:
        raise RuntimeError(f"Cannot find #CHROM header in {vcf_file}.")
    return header_lines, ids, contig_lengths

def prepare_header_context(ids, aids, lists, list_requested, keep_ancestral):
    list_seen = set(lists)
    aid_seen = set(aids)
    aid_indices = []
    for i in range(9, len(ids)):
        if ids[i] in aid_seen:
            aid_indices.append(i)
    if not aid_indices:
        raise RuntimeError("Cannot find any ancestral ID in the VCF.")

    aid_index_seen = set(aid_indices)
    sample_indices = []
    sample_names = []
    for k in range(9, len(ids)):
        sample = ids[k]
        if not keep_ancestral and k in aid_index_seen:
            continue
        if list_requested and sample not in list_seen:
            continue
        sample_indices.append(k)
        sample_names.append(sample)
    fixed_fields = ids[:9]
    return aid_indices, sample_indices, sample_names, fixed_fields

def get_tabix_contigs(vcf_file):
    try:
        contigs = subprocess.run(["tabix", "-l", vcf_file], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Cannot find contigs from tabix index for {vcf_file}: {e}") from e
    return contigs.stdout.strip().splitlines()

def get_ancestral_allele(elements, aid_indices):
    info = elements[8].split(":")
    nucls = elements[4].split(",")
    nucls.insert(0, elements[3])
    anc_alleles = []
    for aid in aid_indices:
        if str(aid).isdigit():
            anc_alleles.append(elements[int(aid)])
    ad = None
    for k, field in enumerate(info):
        if field == "AD":
            ad = k
            break
    anc_nucl_array = []
    for anc_allele in anc_alleles:
        anc_info = anc_allele.split(":")
        alleles = re.split(r"[/|]", anc_info[0])
        if ad is not None and ad < len(anc_info):
            depths = anc_info[ad].split(",")
            depths = [0 if d == "." else int(d) for d in depths]
            if len(alleles) == 2:
                sorted_depths = sorted(
                    range(len(depths)),
                    key=lambda m: depths[m],
                    reverse=True,
                )
                if depths[sorted_depths[0]] == 0:
                    anc_number = "N"
                else:
                    sorted_alleles = [alleles[m] for m in sorted_depths]
                    anc_number = sorted_alleles[0]
            else:
                anc_number = alleles[0]

            if anc_number == "N":
                anc_nucl = anc_number
            else:
                anc_nucl = nucls[int(anc_number)]
        else:
            anc_number = alleles[0]
            anc_nucl = nucls[int(anc_number)]
        anc_nucl_array.append(anc_nucl)
    counts = {}
    first_seen = {}
    for i, nucl in enumerate(anc_nucl_array):
        counts[nucl] = counts.get(nucl, 0) + 1
        if nucl not in first_seen:
            first_seen[nucl] = i
    return min(counts, key=lambda k: (-counts[k], first_seen[k]))

def rename_chr(name, rename_list):
    ori_name = name
    for list in rename_list:
        elements = re.split(r"\s+", list.strip())
        if elements[0] == name:
            name = elements[1]
            break
    if name == ori_name:
        print(f"{ori_name} is unchanged.")
    return name

def process_recode_line(
    line,
    aid_indices,
    sample_indices,
    rename_list,
    biallele,
    hap,
    no_missing,
    shapeit_format,
    thap,):
    elements = line.rstrip("\r\n").split("\t")
    nucl = get_ancestral_allele(elements, aid_indices).upper()
    elements[2] = f"{elements[0]}_{elements[1]}"
    check_existance = f"{elements[3]},{elements[4]}"
    if "*" in nucl:
        if "*" not in check_existance:
            return None, None
        nucl = nucl.replace("*", "B")
    if "*" in check_existance:
        check_existance = check_existance.replace("*", "B")
    check_alleles = check_existance.split(",")
    if nucl not in check_alleles:
        return None, None
    nucl = nucl.replace("B", "*")
    check_existance = check_existance.replace("B", "*")
    sorted_nucls = check_existance.split(",")
    anc = None
    derivs = []
    for i, sorted_nucl in enumerate(sorted_nucls):
        if nucl == sorted_nucl:
            anc = i
        else:
            derivs.append(i)
    if anc is None:
        return None, None
    derivs.insert(0, anc)
    sorted_nucls = [sorted_nucls[i] for i in derivs]
    if biallele and not hap and len(derivs) > 2:
        return None, None
    out_elements = []
    for sample_index in sample_indices:
        out_alleles = []
        gt = elements[sample_index].split(":")
        gt_sep = "|" if "|" in gt[0] else "/"
        alleles = re.split(r"[/|]", gt[0])
        if len(alleles) != 2 and not hap:
            if not no_missing:
                alleles = [".", "."]
            else:
                if len(alleles) == 1:
                    alleles = [alleles[0], alleles[0]]
                else:
                    alleles = ["0", "0"]
            out_elements.append(gt_sep.join(alleles))
        elif len(alleles) == 1 and hap:
            if alleles[0].isdigit():
                for k, deriv in enumerate(derivs):
                    if int(alleles[0]) == deriv:
                        out_alleles.append(str(k))
                        break
                out_elements.append(out_alleles[0] if out_alleles else ".")
            else:
                out_elements.append("0" if no_missing else ".")
        else:
            for l in range(2):
                if alleles[l].isdigit():
                    for k, deriv in enumerate(derivs):
                        if int(alleles[l]) == deriv:
                            out_alleles.append(str(k))
                            break
                else:
                    out_alleles.append("0" if no_missing else ".")
            if not hap:
                out_elements.append(gt_sep.join(out_alleles))
            else:
                out_elements.append(out_alleles[0])
    alt_nucls = sorted_nucls[1:]
    for i in range(9):
        if i == 0 and rename_list:
            elements[i] = rename_chr(elements[i], rename_list)
        if i == 3:
            elements[i] = "0" if thap else sorted_nucls[0]
        if i == 4:
            if not biallele:
                if thap:
                    elements[i] = ",".join(str(j + 1) for j in range(len(alt_nucls)))
                else:
                    elements[i] = ",".join(alt_nucls)
            else:
                if thap:
                    elements[i] = "1"
        if i == 8:
            elements[i] = "GT"
    if biallele:
        unique = sorted(
            {
                allele
                for gt in out_elements
                for allele in re.split(r"[/|]", str(gt))
                if allele != "."
            },
            key=int,
        )
        if len(unique) != 2:
            return None, None
        if not thap:
            elements[4] = alt_nucls[int(unique[1]) - 1]
    map_line = None
    if shapeit_format or thap:
        if shapeit_format:
            out_line = f"{elements[0]} {elements[2]} {elements[1]} {elements[3]} {elements[4]} "
        else:
            out_line = ""
            map_line = f"{elements[2]} {elements[0]} {elements[1]} {elements[3]} {elements[4]}"
        joint_out = " ".join(map(str, out_elements))
        joint_out = re.sub(r"[/|]", " ", joint_out)
        out_line = f"{out_line}{joint_out}"
    else:
        fixed = elements[:9]
        out_line = "\t".join(map(str, fixed + out_elements))
    return out_line, map_line


def process_region(vcf,chr,region_start,region_end,chunk_file,map_chunk_file,aid_indices,sample_indices,rename_list,biallele,hap,no_missing,shapeit_format,thap):
    region = f"{chr}:{region_start}-{region_end}" if region_end else f"{chr}:{region_start}-"
    written_lines = 0
    try:
        with subprocess.Popen(
                ["tabix", "-h", vcf, region],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
        ) as proc:
            assert proc.stdout is not None
            with ExitStack() as stack:
                chunk_out = stack.enter_context(
                    io.TextIOWrapper(pysam.BGZFile(chunk_file, "w"), encoding="utf-8")
                )
                map_out = None
                if thap:
                    map_out = stack.enter_context(
                        io.TextIOWrapper(pysam.BGZFile(map_chunk_file, "w"), encoding="utf-8")
                    )
                for line in proc.stdout:
                    if line.startswith("#"):
                        continue
                    line = line.rstrip("\r\n")
                    out_line, map_line = process_recode_line(line,aid_indices,sample_indices,rename_list,biallele,hap,no_missing,shapeit_format,thap)
                    if out_line is not None:
                        chunk_out.write(out_line + "\n")
                        written_lines += 1
                    if map_out is not None and map_line is not None:
                        map_out.write(map_line + "\n")
            stderr = proc.stderr.read() if proc.stderr is not None else ""
            retcode = proc.wait()
            if retcode != 0:
                raise RuntimeError(f"Cannot process region {region}: {stderr}")
        return written_lines
    except Exception as e:
        raise RuntimeError(f"Cannot process region {region}: {e}") from e

def stream_gzip_file(path, out_handle):
    with pysam.BGZFile(path, "r") as raw_in:
        with io.TextIOWrapper(raw_in, encoding="utf-8") as fh:
            for line in fh:
                out_handle.write(line)

def run_shapeit(unphased_vcf,tagged_vcf,phased_bcf,genetic_map,threads,shapeit_conda):
    try:
        subprocess.run(["tabix", "-f", "-p", "vcf", unphased_vcf], check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"tabix failed with exit code {exc.returncode}.") from exc
    try:
        subprocess.run(["bcftools", "+fill-tags", unphased_vcf, "-Oz", "-o", tagged_vcf, "--", "-t", "AC,AN"], check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Cannot add AC/AN tags to {unphased_vcf} with bcftools +fill-tags.") from exc
    try:
        subprocess.run(["tabix", "-f", "-p", "vcf", tagged_vcf], check=True,)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Cannot index {tagged_vcf} with tabix.") from exc
    contigs = get_tabix_contigs(tagged_vcf)
    if len(contigs) != 1:
        raise RuntimeError(f"SHAPEIT5 phasing expects one chromosome per converter run; found: {','.join(map(str,contigs))}")
    try:
        subprocess.run(["conda", "run", "-n", shapeit_conda, "SHAPEIT5_phase_common", "--input", tagged_vcf, "--map", genetic_map, "--region", contigs[0], "--output", phased_bcf, "--thread", str(threads)], check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"SHAPEIT5_phase_common failed for {tagged_vcf}.") from exc

def convert_phased_bcf_to_haps(phased_bcf, haps_path, sample_path, output_hap):
    try:
        with open(haps_path, "w", encoding="utf-8") as hap_out, \
             open(sample_path, "w", encoding="utf-8") as sample, \
             subprocess.Popen(
                 ["bcftools", "view", phased_bcf],
                 text=True,
                 stdout=subprocess.PIPE,
                 stderr=subprocess.PIPE,
                 encoding="utf-8",
             ) as proc:
            assert proc.stdout is not None
            sample_names = []
            for line in proc.stdout:
                line = line.rstrip("\r\n")
                if line.startswith("#CHROM"):
                    fields = line.split("\t")
                    sample_names = fields[9:]
                    sample.write("ID_1 ID_2 missing\n0 0 0\n")
                    for sample_name in sample_names:
                        if output_hap:
                            sample.write(f"{sample_name} NA 0\n")
                        else:
                            sample.write(f"{sample_name} {sample_name} 0\n")
                    continue

                if line.startswith("#"):
                    continue
                if not sample_names:
                    raise RuntimeError(f"Cannot find sample header in phased BCF {phased_bcf}.")
                elements = line.split("\t")
                alts = elements[4].split(",")
                hap_alleles = []
                for element in elements[9:]:
                    genotype = element.split(":", maxsplit=1)[0]
                    alleles = re.split(r"[/|]", genotype)
                    if len(alleles) == 1:
                        alleles.append(alleles[0])
                    alleles_to_output = [alleles[0]] if output_hap else alleles[:2]
                    for allele in alleles_to_output:
                        if not allele or allele == ".":
                            hap_alleles.append(".")
                        elif allele.isdigit() and int(allele) > 1:
                            hap_alleles.append("1")
                        else:
                            hap_alleles.append(allele)
                output_fields = [elements[0], elements[2], elements[1], elements[3], alts[0], *hap_alleles]
                hap_out.write(" ".join(output_fields) + "\n")
            return_code = proc.wait()
            if return_code != 0:
                stderr = proc.stderr.read().strip() if proc.stderr else ""
                message = f"bcftools view failed for phased BCF {phased_bcf}."
                if stderr:
                    message += f" {stderr}"
                raise RuntimeError(message)
    except FileNotFoundError as exc:
        raise RuntimeError("bcftools was not found in PATH.") from exc

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = "This script convert vcf to an ancestral allele polarized vcf."
    )
    parser.add_argument("-v","--vcf", required=True, help="Input a VCF file.")
    parser.add_argument("-aid","--ancestralID", required=True, help="Input a ancestral ID.")
    parser.add_argument("-k","--keep", action="store_true")
    parser.add_argument("-hp","--hap", action="store_true")
    parser.add_argument("-bi","--biallele", action="store_true")
    parser.add_argument("-nm","--no_missing", action="store_true")
    parser.add_argument("-rchr","--rename_chr")
    parser.add_argument("-l","--list")
    parser.add_argument("-sf","--shapeit_format", action="store_true")
    parser.add_argument("-thap","--thap", action="store_true")
    parser.add_argument("-sc","--shapeit_conda", default="shapeit")
    parser.add_argument("-m","--map")
    parser.add_argument("-o","--out")
    parser.add_argument("-t","--threads", default=1, type=int)
    return parser, parser.parse_args()

def main():
    parser, args = parse_arguments()
    rename_list = []
    aids = []
    missing = []
    TMP_CHUNKS = []
    run_shapeit_requested = args.shapeit_format
    output_hap = args.hap
    if args.shapeit_format:
        if not args.shapeit_conda:
            missing.append("--shapeit_conda")
        if not args.map:
            missing.append("--map")
        if missing:
            parser.error(f"The following arguments are required: {', '.join(missing)}")

    if not args.out:
        args.out = args.vcf.removesuffix(".gz")
        if args.list:
            listname = Path(args.list).name
            listname = listname.removesuffix(".txt").removesuffix(".list")
            args.out = f"{args.out.removesuffix('vcf')}{listname}.vcf"
        args.out = f"{args.out.removesuffix('vcf')}anc.vcf.gz"
        if args.shapeit_format:
            args.out = f"{args.out.removesuffix('vcf.gz')}haps"
        if args.thap:
            args.out = f"{args.out.removesuffix('vcf.gz')}thap"
    else:
        if args.thap:
            args.out = f"{args.out}.thap"
        elif not args.shapeit_format:
            args.out = f"{args.out}.vcf.gz"
        else:
            args.out = f"{args.out}.haps"

    if run_shapeit_requested:
        shapeit_haps_out = args.out
        shapeit_sample_out = f"{args.out.removesuffix('.haps')}.sample"
        shapeit_unphased_vcf = f"{args.out}.unphased.vcf.gz"
        shapeit_tagged_vcf = f"{args.out}.shapeit5_input.vcf.gz"
        shapeit_phased_bcf = f"{args.out}.shapeit5_phased.bcf"
        TMP_CHUNKS.extend([shapeit_unphased_vcf,
                           f"{shapeit_unphased_vcf}.tbi",
                           shapeit_tagged_vcf,
                           f"{shapeit_tagged_vcf}.tbi",
                           shapeit_phased_bcf,
                           f"{shapeit_phased_bcf}.csi"])
        args.out = shapeit_unphased_vcf
        args.shapeit_format = False
        args.hap = False

    map_output = f"{args.out.removesuffix('.thap')}.map" if args.thap else None

    try:
        if args.out.endswith(".gz"):
            out = io.TextIOWrapper(
                pysam.BGZFile(args.out, "w", index=None), encoding="utf-8"
            )
        else:
            out = open(args.out, "w", encoding="utf-8")
        sample = open(map_output, "w", encoding="utf-8") if map_output else None
    except OSError as e:
        parser.error(f"Cannot write output file: {e}")

    print("Start processing VCF file")

    try:
        header_lines, ids, contig_lengths = read_vcf_header(args.vcf)
    except (OSError, RuntimeError) as e:
        parser.error(str(e))

    lists = []
    if args.list:
        try:
            with open(args.list, "r") as f:
                lists = [line.strip() for line in f]
                if args.list.endswith("poplabels"):
                    lists = lists[2:]
                lists = [line.split()[0] for line in lists if line.split()]
        except OSError as e:
            parser.error(f"Cannot read {args.list}: {e}")

    ancestral_path = Path(args.ancestralID)
    if ancestral_path.is_file():
        try:
            with ancestral_path.open("r", encoding="utf-8") as f:
                aids = [line.split()[0] for line in f if line.split()]
        except OSError as e:
            parser.error(f"Cannot read {args.ancestralID}: {e}")
    else:
        aids = [aid.strip() for aid in args.ancestralID.split(",") if aid.strip()]

    try:
        context = prepare_header_context(ids, aids, lists, bool(args.list), args.keep)
    except (OSError, RuntimeError) as e:
        parser.error(str(e))

    aid_indices = context[0]
    sample_indices = context[1]
    sample_names = context[2]
    fixed_fields = context[3]
    if not args.shapeit_format and not args.thap:
        for header_line in header_lines:
            out.write(f"{header_line}\n")
        out.write("\t".join(fixed_fields + sample_names) + "\n")

    if args.rename_chr:
        try:
            with open(args.rename_chr, mode="r", encoding="utf-8") as f:
                rename_list = [line.strip() for line in f]
        except OSError as e:
            parser.error(f"Cannot read {args.rename_chr}: {e}")

    if args.threads > 1:
        if not args.vcf.endswith(".gz"):
            parser.error("-n > 1 requires bgzip-compressed VCF input ending in .gz.")
        if not Path(f"{args.vcf}.tbi").exists():
            try:
                subprocess.run(
                    ["tabix", args.vcf],
                    check=True,
                )
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                parser.error(f"Cannot index {args.vcf} with tabix: {e}")
        try:
            contigs = get_tabix_contigs(args.vcf)
        except (OSError, RuntimeError) as e:
            parser.error(str(e))
        for contig in contigs:
            contig = contig.strip()
            chr_length = contig_lengths.get(contig)
            thread_for_chr = args.threads

            if not chr_length:
                known_contigs = sorted(contig_lengths.keys())
                known = ",".join(known_contigs[:10]) if known_contigs else "none"
                print(
                    f"Cannot find contig length for {contig} in the VCF header."
                    f"Parsed contigs with lengths: {known}."
                    f"Processing this contig as one whole region instead of splitting across threads.",
                    file=sys.stderr,
                )
                thread_for_chr = 1
            else:
                thread_for_chr = min(thread_for_chr, chr_length)

            chunk_files = []
            map_chunk_files = []
            chunk_regions = []
            futures = []
            with ProcessPoolExecutor(max_workers=thread_for_chr) as executor:
                for i in range(1, thread_for_chr + 1):
                    if not chr_length:
                        region_start = 1
                        region_end = 0
                    elif i < thread_for_chr:
                        if i == 1:
                            region_start = 1
                        else:
                            region_start = int(chr_length / thread_for_chr * (i - 1)) + 1
                        region_end = int(chr_length / thread_for_chr * i)
                    else:
                        region_start = int(chr_length / thread_for_chr * (i - 1)) + 1
                        region_end = 0

                    chunk_chr = re.sub(r"[^A-Za-z0-9_.-]", "_", contig)

                    tmp = tempfile.NamedTemporaryFile(
                        prefix=f"vcf2anc.{chunk_chr}.{i}.",
                        suffix=".gz",
                        dir=Path(args.vcf).parent,
                        delete=False,
                    )
                    chunk_file = tmp.name
                    tmp.close()
                    TMP_CHUNKS.append(chunk_file)
                    chunk_files.append(chunk_file)

                    if args.thap:
                        tmp = tempfile.NamedTemporaryFile(
                            prefix=f"vcf2anc.{chunk_chr}.{i}.map.",
                            suffix=".gz",
                            dir=Path(args.vcf).parent,
                            delete=False,
                        )
                        map_file = tmp.name
                        tmp.close()
                        TMP_CHUNKS.append(map_file)
                        map_chunk_files.append(map_file)
                    else:
                        map_chunk_files.append(None)

                    chunk_regions.append(f"{contig}:{region_start}-{region_end}" if region_end else f"{contig}:{region_start}-")

                    futures.append(
                        executor.submit(
                            process_region,
                            args.vcf,
                            contig,
                            region_start,
                            region_end,
                            chunk_file,
                            map_chunk_files[-1],
                            aid_indices,
                            sample_indices,
                            rename_list,
                            args.biallele,
                            args.hap,
                            args.no_missing,
                            args.shapeit_format,
                            args.thap,
                        )
                    )

                for future in futures:
                    future.result()

            if out is None:
                raise RuntimeError("Output handle was not initialized.")

            for idx, chunk_file in enumerate(chunk_files):
                stream_gzip_file(chunk_file, out)
                try:
                    Path(chunk_file).unlink()
                except FileNotFoundError:
                    pass
                except OSError as e:
                    print(f"Cannot remove temporary output {chunk_file}: {e}", file=sys.stderr)
                try:
                    TMP_CHUNKS.remove(chunk_file)
                except ValueError:
                    pass

                map_file = map_chunk_files[idx]
                if args.thap and sample is not None and map_file is not None:
                    stream_gzip_file(map_file, sample)
                    try:
                        Path(map_file).unlink()
                    except FileNotFoundError:
                        pass
                    except OSError as e:
                        print(f"Cannot remove temporary output {map_file}: {e}", file=sys.stderr)
                    try:
                        TMP_CHUNKS.remove(map_file)
                    except ValueError:
                        pass
    else:
        if args.vcf.endswith(".gz"):
            binary_in = pysam.BGZFile(args.vcf, "r")
            raw_in = io.TextIOWrapper(binary_in)
        else:
            raw_in = open(args.vcf, "r")
        with raw_in:
            for line in raw_in:
                if line.startswith("#"):
                    continue
                line = line.rstrip("\r\n")
                out_line, map_line = process_recode_line(line,aid_indices,sample_indices,rename_list,args.biallele,args.hap,args.no_missing,args.shapeit_format,args.thap)
                if not out_line:
                    continue
                out.write(f'{out_line}\n')
                if args.thap and sample is not None and map_line is not None:
                    sample.write(f'{map_line}\n')
    out.close()
    if sample is not None:
        sample.close()

    #run shapeit5
    if run_shapeit_requested:
        run_shapeit(shapeit_unphased_vcf,shapeit_tagged_vcf,shapeit_phased_bcf,args.map,args.threads,args.shapeit_conda)
        convert_phased_bcf_to_haps(
            shapeit_phased_bcf,
            shapeit_haps_out,
            shapeit_sample_out,
            output_hap,
        )

    for temporary_file in TMP_CHUNKS:
        try:
            Path(temporary_file).unlink()
        except FileNotFoundError:
            pass
        except OSError as e:
            print(f"Cannot remove temporary file {temporary_file}: {e}", file=sys.stderr)

    print("Done.")


if __name__ == "__main__":
    main()
