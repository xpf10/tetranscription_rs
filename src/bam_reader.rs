/// BAM/SAM parsing and core counting loop using noodles (pure Rust).
/// Ported from bin/TEcount. Uses Rayon for parallel overlap queries.

use std::collections::HashMap;

use pyo3::prelude::*;
use rayon::prelude::*;

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

// ---------------------------------------------------------------------------
// Helper types and functions
// ---------------------------------------------------------------------------

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
}

struct ParsedAlignment {
    clean_name: String,
    info: ReadInfo,
}

fn parse_record_to_alignment(record: &noodles_bam::Record) -> Option<ParsedAlignment> {
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
    };

    Some(ParsedAlignment { clean_name, info })
}

// ---------------------------------------------------------------------------
// Read group: all alignments sharing the same read name
// ---------------------------------------------------------------------------

/// A read group ready for parallel overlap annotation.
enum ReadGroup {
    Paired {
        read1s: Vec<ReadInfo>,
        read2s: Vec<ReadInfo>,
    },
    Single {
        reads: Vec<ReadInfo>,
    },
}

/// Result of overlap annotation for one read group (no mutable state needed).
struct AnnotResult {
    annot_gene: Vec<Vec<String>>,
    annot_te: Vec<Vec<usize>>,
    is_multi: bool,
    /// For paired unique reads: fragment length info
    frag_len: Option<i64>,
}

fn overlap_annotation_readgroup(
    group: &ReadGroup,
    references: &[String],
    gene_index: &GeneIndex,
    te_index: &TEIndex,
    stranded: &str,
) -> AnnotResult {
    let reads: Vec<(Option<ReadInfo>, Option<ReadInfo>)> = match group {
        ReadGroup::Paired { read1s, read2s } => {
            let is_multi = read1s.len() > 1 || read2s.len() > 1;
            if is_multi {
                let mut pairs = Vec::new();
                if read2s.is_empty() {
                    for r in read1s { pairs.push((Some(r.clone()), None)); }
                } else if read1s.is_empty() {
                    for r in read2s { pairs.push((None, Some(r.clone()))); }
                } else if read2s.len() == read1s.len() {
                    for j in 0..read1s.len() {
                        pairs.push((Some(read1s[j].clone()), Some(read2s[j].clone())));
                    }
                }
                pairs
            } else {
                let r1 = if read1s.len() == 1 { Some(read1s[0].clone()) } else { None };
                let r2 = if read2s.len() == 1 { Some(read2s[0].clone()) } else { None };
                vec![(r1, r2)]
            }
        }
        ReadGroup::Single { reads } => {
            reads.iter().map(|r| (Some(r.clone()), None)).collect()
        }
    };

    let is_multi = reads.len() > 1;
    let mut frag_len = None;
    if !is_multi && !reads.is_empty() {
        if let (Some(ref r1), Some(ref r2)) = (&reads[0].0, &reads[0].1) {
            if r1.is_proper_pair {
                let pos1 = r1.reference_start;
                let pos2 = r2.reference_start;
                frag_len = Some((pos1 - pos2).abs() + r2.query_length);
            }
        }
    }

    let mut annot_gene: Vec<Vec<String>> = Vec::new();
    let mut annot_te: Vec<Vec<usize>> = Vec::new();

    for (r1, r2) in &reads {
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

    AnnotResult { annot_gene, annot_te, is_multi, frag_len }
}

// ---------------------------------------------------------------------------
// Main counting function
// ---------------------------------------------------------------------------

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
    // Phase 1: Read all records from BAM
    let mut reader = noodles_bam::io::reader::Builder::default()
        .build_from_path(bam_path)
        .unwrap_or_else(|e| {
            eprintln!("Error opening BAM file {}: {}", bam_path, e);
            std::process::exit(1);
        });

    let header = reader.read_header().expect("Error reading BAM header");
    let references: Vec<String> = header.reference_sequences()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    eprintln!("Reading BAM file...");
    let mut alignments: Vec<ParsedAlignment> = Vec::new();
    let mut record = noodles_bam::Record::default();
    let mut paired = false;
    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {}
            Err(_) => continue,
        }
        if let Some(aln) = parse_record_to_alignment(&record) {
            if aln.info.is_paired { paired = true; }
            alignments.push(aln);
        }
    }
    eprintln!("Read {} alignments.", alignments.len());

    // Phase 2: Sort by name (parallel if sortByPos, already sorted otherwise)
    if sort_by_pos {
        eprintln!("Sorting by read name (parallel)...");
        alignments.par_sort_by(|a, b| a.clean_name.cmp(&b.clean_name));
    }

    // Phase 3: Group by read name into ReadGroups
    eprintln!("Grouping reads by name...");
    let mut groups: Vec<ReadGroup> = Vec::new();

    if paired {
        let mut i = 0;
        while i < alignments.len() {
            let name = &alignments[i].clean_name;
            let start = i;
            while i < alignments.len() && &alignments[i].clean_name == name {
                i += 1;
            }
            let mut read1s = Vec::new();
            let mut read2s = Vec::new();
            for j in start..i {
                let info = alignments[j].info.clone();
                if info.is_read1 { read1s.push(info.clone()); }
                if info.is_read2 { read2s.push(info); }
            }
            groups.push(ReadGroup::Paired { read1s, read2s });
        }
    } else {
        let mut i = 0;
        while i < alignments.len() {
            let name = &alignments[i].clean_name;
            let start = i;
            while i < alignments.len() && &alignments[i].clean_name == name {
                i += 1;
            }
            let reads: Vec<ReadInfo> = (start..i).map(|j| alignments[j].info.clone()).collect();
            groups.push(ReadGroup::Single { reads });
        }
    }
    eprintln!("{} read groups formed.", groups.len());

    // Phase 4: Parallel overlap annotation (read-only on indexes)
    eprintln!("Running parallel overlap annotation...");
    let annot_results: Vec<AnnotResult> = groups
        .par_iter()
        .map(|group| overlap_annotation_readgroup(group, &references, gene_index, te_index, stranded))
        .collect();

    // Phase 5: Serial count aggregation
    eprintln!("Aggregating counts...");
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

    for res in &annot_results {
        if res.is_multi {
            nonunique += 1;
            if te_mode == "uniq" {
                empty += 1;
                continue;
            }
        } else {
            uniq_reads += 1;
            if let Some(fl) = res.frag_len {
                if fl <= max_length && tmp_cnt < 10000 {
                    avg_read_length += fl;
                    tmp_cnt += 1;
                }
            }
        }

        let num_alignments = if res.is_multi {
            res.annot_gene.len().max(res.annot_te.len()).max(1)
        } else {
            1
        };

        if num_alignments > 1 {
            let no_annot_te = parse_annotations_te(
                &res.annot_te, &mut te_counts, &mut te_multi_counts,
                &mut multi_reads, &mut leftover_te,
            );
            if no_annot_te {
                let no_annot_gene = parse_annotations_gene(
                    &res.annot_gene, &mut gene_counts, &mut leftover_gene,
                );
                if no_annot_gene { empty += 1; }
            }
        } else {
            let no_annot_gene = parse_annotations_gene(
                &res.annot_gene, &mut gene_counts, &mut leftover_gene,
            );
            if no_annot_gene {
                let no_annot_te = parse_annotations_te(
                    &res.annot_te, &mut te_counts, &mut te_multi_counts,
                    &mut multi_reads, &mut leftover_te,
                );
                if no_annot_te { empty += 1; }
            }
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
