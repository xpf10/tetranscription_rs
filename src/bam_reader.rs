/// BAM/SAM parsing and core counting loop using noodles (pure Rust).
/// Ported from bin/TEcount.

use std::collections::HashMap;

use pyo3::prelude::*;

use crate::annotation::{
    parse_annotations_gene, parse_annotations_te, resolve_annotation_ambiguity_gene,
    resolve_annotation_ambiguity_te,
};
use crate::em_algorithm::em_estimate;
use crate::gene_index::GeneIndex;
use crate::te_index::TEIndex;
use crate::types::{CigarElement, ExonInterval, Strand};

/// Python-facing count result
#[pyclass]
pub struct PyCountResult {
    #[pyo3(get)]
    pub gene_counts: HashMap<String, f64>,
    #[pyo3(get)]
    pub te_instance_counts: Vec<f64>,
    #[pyo3(get)]
    pub te_element_counts: HashMap<String, f64>,
    #[pyo3(get)]
    pub total_annotated: i64,
    #[pyo3(get)]
    pub total_nonunique: i64,
    #[pyo3(get)]
    pub total_unannotated: i64,
}

#[pyfunction(name = "count_transcript_abundance")]
pub fn count_transcript_abundance_py(
    bam_path: &str,
    gene_index: &GeneIndex,
    te_index: &TEIndex,
    stranded: &str,
    te_mode: &str,
    sort_by_pos: bool,
    num_iterations: i32,
    frag_length: i64,
    max_length: i64,
) -> PyResult<PyCountResult> {
    let result = count_transcript_abundance(
        bam_path, gene_index, te_index, stranded, te_mode,
        sort_by_pos, num_iterations, frag_length, max_length,
    );
    let te_element_counts = te_index.group_by_element(result.te_instance_counts.clone());
    Ok(PyCountResult {
        gene_counts: result.gene_counts,
        te_instance_counts: result.te_instance_counts,
        te_element_counts,
        total_annotated: result.total_annotated,
        total_nonunique: result.total_nonunique,
        total_unannotated: result.total_unannotated,
    })
}

struct CountResult {
    gene_counts: HashMap<String, f64>,
    te_instance_counts: Vec<f64>,
    total_annotated: i64,
    total_nonunique: i64,
    total_unannotated: i64,
}

fn fetch_exon(chrom: &str, pos: i64, cigar: &[CigarElement], direction: i32) -> Vec<ExonInterval> {
    let mut chrom_st = pos + 1;
    let mut result = Vec::new();
    let strand = match direction {
        1 => Strand::Plus,
        -1 => Strand::Minus,
        _ => Strand::Unknown,
    };
    for c in cigar {
        match c.code {
            0 => {
                result.push(ExonInterval {
                    chrom: chrom.to_string(),
                    start: chrom_st,
                    end: chrom_st + c.len - 1,
                    strand: strand.clone(),
                });
                chrom_st += c.len;
            }
            2 | 3 | 4 => { chrom_st += c.len; }
            _ => {}
        }
    }
    result
}

fn get_direction(r1_reverse: Option<bool>, r2_reverse: Option<bool>, stranded: &str) -> i32 {
    let mut direction = 1;
    if r1_reverse == Some(true) { direction = -1; }
    if r2_reverse == Some(false) { direction = -1; }
    match stranded {
        "no" => 0,
        "reverse" => direction * -1,
        _ => direction,
    }
}

fn convert_cigar(record: &noodles_bam::Record) -> Vec<CigarElement> {
    let mut result = Vec::new();
    for op_result in record.cigar().iter() {
        let op_val = match op_result {
            Ok(v) => v,
            Err(_) => continue,
        };
        let len: i64 = op_val.len() as i64;
        let code: u32 = if op_val.kind() == noodles_sam::alignment::record::cigar::op::Kind::Match { 0 }
            else if op_val.kind() == noodles_sam::alignment::record::cigar::op::Kind::Insertion { 1 }
            else if op_val.kind() == noodles_sam::alignment::record::cigar::op::Kind::Deletion { 2 }
            else if op_val.kind() == noodles_sam::alignment::record::cigar::op::Kind::Skip { 3 }
            else if op_val.kind() == noodles_sam::alignment::record::cigar::op::Kind::SoftClip { 4 }
            else if op_val.kind() == noodles_sam::alignment::record::cigar::op::Kind::HardClip { 5 }
            else if op_val.kind() == noodles_sam::alignment::record::cigar::op::Kind::SequenceMatch { 7 }
            else if op_val.kind() == noodles_sam::alignment::record::cigar::op::Kind::SequenceMismatch { 8 }
            else { continue; };
        result.push(CigarElement { code, len });
    }
    result
}

#[derive(Clone)]
struct ReadInfo {
    tid: usize,
    pos: i64,
    cigar: Vec<CigarElement>,
    is_reverse: bool,
    is_read1: bool,
    is_read2: bool,
    is_paired: bool,
    is_proper_pair: bool,
    query_length: i64,
    reference_start: i64,
    clean_name: String,
}

/// A filtered, parsed alignment ready for grouping by name.
struct ParsedAlignment {
    clean_name: String,
    info: ReadInfo,
}

fn parse_record_to_alignment(record: &noodles_bam::Record, _references: &[String]) -> Option<ParsedAlignment> {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_duplicate() { return None; }
    if flags.bits() & 512 != 0 { return None; }

    let cur_read_name = record.name().map(|n| n.to_string()).unwrap_or_default();
    let cigar = convert_cigar(record);
    let tid = match record.reference_sequence_id() {
        Some(Ok(id)) => id,
        _ => return None,
    };
    let pos = match record.alignment_start() {
        Some(Ok(p)) => p.get() as i64 - 1,
        _ => return None,
    };
    let seq_len = record.sequence().len() as i64;
    let is_paired = flags.is_segmented();

    let clean_name = if is_paired {
        let added_flag_pos = cur_read_name.find('/').unwrap_or(cur_read_name.len());
        cur_read_name[..added_flag_pos].to_string()
    } else {
        cur_read_name
    };

    let info = ReadInfo {
        tid,
        pos,
        cigar,
        is_reverse: flags.is_reverse_complemented(),
        is_read1: flags.is_first_segment(),
        is_read2: flags.is_last_segment(),
        is_paired,
        is_proper_pair: flags.is_properly_segmented(),
        query_length: seq_len,
        reference_start: pos,
        clean_name: clean_name.clone(),
    };

    Some(ParsedAlignment { clean_name, info })
}

fn count_transcript_abundance(
    bam_path: &str,
    gene_index: &GeneIndex,
    te_index: &TEIndex,
    stranded: &str,
    te_mode: &str,
    sort_by_pos: bool,
    num_iterations: i32,
    frag_length: i64,
    max_length: i64,
) -> CountResult {
    let mut reader = noodles_bam::io::reader::Builder::default()
        .build_from_path(bam_path)
        .unwrap_or_else(|e| {
            eprintln!("Error opening BAM file {}: {}", bam_path, e);
            std::process::exit(1);
        });

    let header = reader.read_header()
        .expect("Error reading BAM header");

    let references: Vec<String> = header.reference_sequences()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    let mut gene_counts: HashMap<String, f64> = HashMap::new();
    for f in gene_index.features() {
        gene_counts.insert(f, 0.0);
    }

    let n_te = te_index.num_instances();
    let mut te_counts = vec![0.0; n_te];
    let mut te_multi_counts = vec![0.0; n_te];
    let mut multi_reads: Vec<Vec<usize>> = Vec::new();
    let mut leftover_gene: Vec<(Vec<Vec<String>>, f64)> = Vec::new();
    let mut leftover_te: Vec<(Vec<Vec<usize>>, f64)> = Vec::new();

    let mut empty: i64 = 0;
    let mut nonunique: i64 = 0;
    let mut uniq_reads: i64 = 0;
    let mut avg_read_length: i64 = 0;
    let mut tmp_cnt: i64 = 0;
    let mut paired = false;

    if sort_by_pos {
        // ---- sortByPos path: collect all records, sort by name, then group ----
        eprintln!("BAM sorted by position. Reading all alignments into memory for name-sorting...");

        let mut alignments: Vec<ParsedAlignment> = Vec::new();
        let mut record = noodles_bam::Record::default();
        loop {
            match reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {}
                Err(_) => continue,
            }
            if let Some(aln) = parse_record_to_alignment(&record, &references) {
                if aln.info.is_paired { paired = true; }
                alignments.push(aln);
            }
        }
        eprintln!("Read {} alignments. Sorting by read name...", alignments.len());
        alignments.sort_by(|a, b| a.clean_name.cmp(&b.clean_name));
        eprintln!("Sorting done. Processing...");

        let mut i: i64 = 0;
        let mut prev_name = String::new();
        let mut multi_read1: Vec<ReadInfo> = Vec::new();
        let mut multi_read2: Vec<ReadInfo> = Vec::new();
        let mut alignments_per_read: Vec<(Option<ReadInfo>, Option<ReadInfo>)> = Vec::new();

        for aln in alignments {
            i += 1;
            let clean_name = aln.clean_name;
            let read_info = aln.info;

            if paired {
                if clean_name == prev_name || prev_name.is_empty() {
                    prev_name = clean_name;
                    if read_info.is_read1 { multi_read1.push(read_info.clone()); }
                    if read_info.is_read2 { multi_read2.push(read_info); }
                    continue;
                }

                process_read_group_paired(
                    &multi_read1, &multi_read2, &references, gene_index, te_index,
                    stranded, te_mode, &mut alignments_per_read,
                    &mut gene_counts, &mut te_counts, &mut te_multi_counts,
                    &mut multi_reads, &mut leftover_gene, &mut leftover_te,
                    &mut uniq_reads, &mut nonunique, &mut empty,
                    &mut avg_read_length, &mut tmp_cnt, max_length,
                );
                alignments_per_read.clear();
                multi_read1.clear();
                multi_read2.clear();
                prev_name = clean_name;
                if read_info.is_read1 { multi_read1.push(read_info.clone()); }
                if read_info.is_read2 { multi_read2.push(read_info); }
            } else {
                if clean_name == prev_name || prev_name.is_empty() {
                    alignments_per_read.push((Some(read_info), None));
                    prev_name = clean_name;
                    continue;
                }

                if tmp_cnt < 10000 {
                    if let Some(ref r) = alignments_per_read[0].0 {
                        avg_read_length += r.query_length;
                        tmp_cnt += 1;
                    }
                }

                if alignments_per_read.len() == 1 {
                    uniq_reads += 1;
                } else {
                    nonunique += 1;
                    if te_mode == "uniq" {
                        empty += 1;
                        alignments_per_read.clear();
                        prev_name = clean_name;
                        alignments_per_read.push((Some(read_info), None));
                        continue;
                    }
                }

                annotate_and_count(
                    &alignments_per_read, &references, gene_index, te_index,
                    stranded, &mut gene_counts, &mut te_counts, &mut te_multi_counts,
                    &mut multi_reads, &mut leftover_gene, &mut leftover_te, &mut empty,
                );

                if i % 1_000_000 == 0 {
                    eprintln!("{} alignments processed.", i);
                }
                alignments_per_read.clear();
                prev_name = clean_name;
                alignments_per_read.push((Some(read_info), None));
            }
        }

        // Process last group
        if paired {
            process_read_group_paired(
                &multi_read1, &multi_read2, &references, gene_index, te_index,
                stranded, te_mode, &mut alignments_per_read,
                &mut gene_counts, &mut te_counts, &mut te_multi_counts,
                &mut multi_reads, &mut leftover_gene, &mut leftover_te,
                &mut uniq_reads, &mut nonunique, &mut empty,
                &mut avg_read_length, &mut tmp_cnt, max_length,
            );
        } else if !alignments_per_read.is_empty() {
            if alignments_per_read.len() == 1 { uniq_reads += 1; } else { nonunique += 1; }
            annotate_and_count(
                &alignments_per_read, &references, gene_index, te_index,
                stranded, &mut gene_counts, &mut te_counts, &mut te_multi_counts,
                &mut multi_reads, &mut leftover_gene, &mut leftover_te, &mut empty,
            );
        }
    } else {
        // ---- Name-sorted BAM path: streaming ----
        let mut i: i64 = 0;
        let mut prev_read_name = String::new();
        let mut alignments_per_read: Vec<(Option<ReadInfo>, Option<ReadInfo>)> = Vec::new();
        let mut multi_read1: Vec<ReadInfo> = Vec::new();
        let mut multi_read2: Vec<ReadInfo> = Vec::new();
        let mut record = noodles_bam::Record::default();

        loop {
            match reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {}
                Err(_) => continue,
            }

            i += 1;
            let Some(aln) = parse_record_to_alignment(&record, &references) else { continue };
            let clean_name = aln.clean_name;
            let read_info = aln.info;

            if read_info.is_paired { paired = true; }

            if paired {
                if clean_name == prev_read_name || prev_read_name.is_empty() {
                    prev_read_name = clean_name;
                    if read_info.is_read1 { multi_read1.push(read_info.clone()); }
                    if read_info.is_read2 { multi_read2.push(read_info); }
                    continue;
                }

                process_read_group_paired(
                    &multi_read1, &multi_read2, &references, gene_index, te_index,
                    stranded, te_mode, &mut alignments_per_read,
                    &mut gene_counts, &mut te_counts, &mut te_multi_counts,
                    &mut multi_reads, &mut leftover_gene, &mut leftover_te,
                    &mut uniq_reads, &mut nonunique, &mut empty,
                    &mut avg_read_length, &mut tmp_cnt, max_length,
                );
                alignments_per_read.clear();
                multi_read1.clear();
                multi_read2.clear();
                prev_read_name = clean_name;
                if read_info.is_read1 { multi_read1.push(read_info.clone()); }
                if read_info.is_read2 { multi_read2.push(read_info); }
            } else {
                if clean_name == prev_read_name || prev_read_name.is_empty() {
                    alignments_per_read.push((Some(read_info), None));
                    prev_read_name = clean_name;
                    continue;
                }

                if tmp_cnt < 10000 {
                    if let Some(ref r) = alignments_per_read[0].0 {
                        avg_read_length += r.query_length;
                        tmp_cnt += 1;
                    }
                }

                if alignments_per_read.len() == 1 {
                    uniq_reads += 1;
                } else {
                    nonunique += 1;
                    if te_mode == "uniq" {
                        empty += 1;
                        alignments_per_read.clear();
                        prev_read_name = clean_name;
                        alignments_per_read.push((Some(read_info), None));
                        continue;
                    }
                }

                annotate_and_count(
                    &alignments_per_read, &references, gene_index, te_index,
                    stranded, &mut gene_counts, &mut te_counts, &mut te_multi_counts,
                    &mut multi_reads, &mut leftover_gene, &mut leftover_te, &mut empty,
                );

                if i % 1_000_000 == 0 {
                    eprintln!("{} alignments processed.", i);
                }
                alignments_per_read.clear();
                prev_read_name = clean_name;
                alignments_per_read.push((Some(read_info), None));
            }
        }

        // Process last read group
        if paired {
            process_read_group_paired(
                &multi_read1, &multi_read2, &references, gene_index, te_index,
                stranded, te_mode, &mut alignments_per_read,
                &mut gene_counts, &mut te_counts, &mut te_multi_counts,
                &mut multi_reads, &mut leftover_gene, &mut leftover_te,
                &mut uniq_reads, &mut nonunique, &mut empty,
                &mut avg_read_length, &mut tmp_cnt, max_length,
            );
        } else if !alignments_per_read.is_empty() {
            if alignments_per_read.len() == 1 { uniq_reads += 1; } else { nonunique += 1; }
            annotate_and_count(
                &alignments_per_read, &references, gene_index, te_index,
                stranded, &mut gene_counts, &mut te_counts, &mut te_multi_counts,
                &mut multi_reads, &mut leftover_gene, &mut leftover_te, &mut empty,
            );
        }
    }

    // Resolve leftover ambiguities
    if !leftover_gene.is_empty() {
        resolve_annotation_ambiguity_gene(&mut gene_counts, &leftover_gene);
    }
    if !leftover_te.is_empty() {
        resolve_annotation_ambiguity_te(&mut te_counts, &leftover_te);
    }

    eprintln!("uniq te counts = {}", te_counts.iter().sum::<f64>() as i64);

    let estimated_read_length = if !paired && frag_length > 0 {
        frag_length
    } else if avg_read_length > 0 && tmp_cnt > 0 {
        avg_read_length / tmp_cnt
    } else {
        100
    };

    let new_te_multi_counts = if num_iterations > 0 && !multi_reads.is_empty() {
        eprintln!("Starting iterative optimization...");
        em_estimate(te_index, &multi_reads, &te_counts, &te_multi_counts,
                    num_iterations, estimated_read_length)
    } else {
        te_multi_counts
    };

    for j in 0..te_counts.len() {
        te_counts[j] += new_te_multi_counts[j];
    }

    let st: f64 = te_counts.iter().sum();
    let sg: f64 = gene_counts.values().sum();

    eprintln!("TE counts total {}", st);
    eprintln!("Gene counts total {}", sg);
    eprintln!("Total annotated = {}", (st + sg) as i64);
    eprintln!("Total non-unique = {}", nonunique);
    eprintln!("Total unannotated = {}", empty);

    CountResult {
        gene_counts,
        te_instance_counts: te_counts,
        total_annotated: (st + sg) as i64,
        total_nonunique: nonunique,
        total_unannotated: empty,
    }
}

#[allow(clippy::too_many_arguments)]
fn process_read_group_paired(
    multi_read1: &[ReadInfo], multi_read2: &[ReadInfo],
    references: &[String], gene_index: &GeneIndex, te_index: &TEIndex,
    stranded: &str, te_mode: &str,
    alignments_per_read: &mut Vec<(Option<ReadInfo>, Option<ReadInfo>)>,
    gene_counts: &mut HashMap<String, f64>,
    te_counts: &mut Vec<f64>, te_multi_counts: &mut Vec<f64>,
    multi_reads: &mut Vec<Vec<usize>>,
    leftover_gene: &mut Vec<(Vec<Vec<String>>, f64)>,
    leftover_te: &mut Vec<(Vec<Vec<usize>>, f64)>,
    uniq_reads: &mut i64, nonunique: &mut i64, empty: &mut i64,
    avg_read_length: &mut i64, tmp_cnt: &mut i64, max_length: i64,
) {
    if multi_read1.len() <= 1 && multi_read2.len() <= 1 {
        *uniq_reads += 1;
        let read1 = if multi_read1.len() == 1 { Some(multi_read1[0].clone()) } else { None };
        let read2 = if multi_read2.len() == 1 { Some(multi_read2[0].clone()) } else { None };

        if let (Some(ref r1), Some(ref r2)) = (&read1, &read2) {
            if r1.is_proper_pair && *tmp_cnt < 10000 {
                let pos1 = r1.reference_start;
                let pos2 = r2.reference_start;
                if (pos1 - pos2).abs() <= max_length {
                    *avg_read_length += (pos1 - pos2).abs() + r2.query_length;
                    *tmp_cnt += 1;
                }
            }
        }
        alignments_per_read.push((read1, read2));
    } else {
        *nonunique += 1;
        if te_mode == "uniq" { *empty += 1; return; }

        if multi_read2.is_empty() {
            for r in multi_read1 { alignments_per_read.push((Some(r.clone()), None)); }
        } else if multi_read1.is_empty() {
            for r in multi_read2 { alignments_per_read.push((None, Some(r.clone()))); }
        } else if multi_read2.len() == multi_read1.len() {
            for j in 0..multi_read1.len() {
                alignments_per_read.push((Some(multi_read1[j].clone()), Some(multi_read2[j].clone())));
            }
        }
    }

    annotate_and_count(
        alignments_per_read, references, gene_index, te_index,
        stranded, gene_counts, te_counts, te_multi_counts,
        multi_reads, leftover_gene, leftover_te, empty,
    );
}

#[allow(clippy::too_many_arguments)]
fn annotate_and_count(
    alignments: &[(Option<ReadInfo>, Option<ReadInfo>)],
    references: &[String], gene_index: &GeneIndex, te_index: &TEIndex,
    stranded: &str,
    gene_counts: &mut HashMap<String, f64>,
    te_counts: &mut Vec<f64>, te_multi_counts: &mut Vec<f64>,
    multi_reads: &mut Vec<Vec<usize>>,
    leftover_gene: &mut Vec<(Vec<Vec<String>>, f64)>,
    leftover_te: &mut Vec<(Vec<Vec<usize>>, f64)>,
    empty: &mut i64,
) {
    let (annot_gene, annot_te) = overlap_annotation(alignments, references, gene_index, te_index, stranded);

    if alignments.len() > 1 {
        let no_annot_te = parse_annotations_te(&annot_te, te_counts, te_multi_counts, multi_reads, leftover_te);
        if no_annot_te {
            let no_annot_gene = parse_annotations_gene(&annot_gene, gene_counts, leftover_gene);
            if no_annot_gene { *empty += 1; }
        }
    } else {
        let no_annot_gene = parse_annotations_gene(&annot_gene, gene_counts, leftover_gene);
        if no_annot_gene {
            let no_annot_te = parse_annotations_te(&annot_te, te_counts, te_multi_counts, multi_reads, leftover_te);
            if no_annot_te { *empty += 1; }
        }
    }
}

fn overlap_annotation(
    reads: &[(Option<ReadInfo>, Option<ReadInfo>)],
    references: &[String], gene_index: &GeneIndex, te_index: &TEIndex,
    stranded: &str,
) -> (Vec<Vec<String>>, Vec<Vec<usize>>) {
    let mut annot_gene: Vec<Vec<String>> = Vec::new();
    let mut annot_te: Vec<Vec<usize>> = Vec::new();

    for (r1, r2) in reads {
        let r1_reverse = r1.as_ref().map(|r| r.is_reverse);
        let r2_reverse = r2.as_ref().map(|r| r.is_reverse);
        let direction = get_direction(r1_reverse, r2_reverse, stranded);

        let mut itv_list: Vec<ExonInterval> = Vec::new();

        if let Some(ref r) = r1 {
            if r.tid < references.len() {
                itv_list.extend(fetch_exon(&references[r.tid], r.pos, &r.cigar, direction));
            }
        }
        if let Some(ref r) = r2 {
            if r.tid < references.len() {
                itv_list.extend(fetch_exon(&references[r.tid], r.pos, &r.cigar, direction));
            }
        }

        let tes = te_index.te_annotation(&itv_list);
        let genes = gene_index.gene_annotation(&itv_list);

        if !tes.is_empty() { annot_te.push(tes); }
        if !genes.is_empty() {
            let mut ug = genes;
            ug.sort();
            ug.dedup();
            annot_gene.push(ug);
        }
    }

    (annot_gene, annot_te)
}
