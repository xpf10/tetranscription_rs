/// BAM/SAM parsing and core counting loop using noodles (pure Rust).
/// Uses batched streaming to bound memory, with Rayon parallelism per batch.

use std::collections::HashMap;
use std::io::Read;

use pyo3::prelude::*;
use rayon::prelude::*;
use smallvec::SmallVec;

use crate::annotation::{
    parse_annotations_gene, parse_annotations_te, resolve_annotation_ambiguity_gene,
    resolve_annotation_ambiguity_te,
};
use crate::em_algorithm::em_estimate;
use crate::gene_index::GeneIndex;
use crate::te_index::TEIndex;
use crate::types::{CigarElement, ExonInterval, Strand};

// ---------------------------------------------------------------------------
// Python-facing types
// ---------------------------------------------------------------------------

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
// Core types
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct ReadInfo {
    tid: usize,
    pos: i64,
    cigar: SmallVec<[CigarElement; 4]>,
    is_reverse: bool,
    is_read1: bool,
    is_read2: bool,
    is_paired: bool,
    is_proper_pair: bool,
    query_length: i64,
    reference_start: i64,
}

enum ReadGroup {
    Paired {
        read1s: Vec<ReadInfo>,
        read2s: Vec<ReadInfo>,
    },
    Single {
        reads: Vec<ReadInfo>,
    },
}

struct AnnotResult {
    annot_gene: Vec<Vec<String>>,
    annot_te: Vec<Vec<usize>>,
    is_multi: bool,
    frag_len: Option<i64>,
}

// ---------------------------------------------------------------------------
// BAM record parsing helpers
// ---------------------------------------------------------------------------

fn convert_cigar(record: &noodles_bam::Record) -> SmallVec<[CigarElement; 4]> {
    let mut result = SmallVec::new();
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

/// Parse a BAM record into (clean_name, ReadInfo). Returns None for unmapped/duplicate/QC-fail.
fn parse_record(record: &noodles_bam::Record) -> Option<(String, ReadInfo)> {
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
        let idx = cur_read_name.find('/').unwrap_or(cur_read_name.len());
        cur_read_name[..idx].to_string()
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

    Some((clean_name, info))
}

// ---------------------------------------------------------------------------
// Exon interval extraction
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

// ---------------------------------------------------------------------------
// Overlap annotation for a read group
// ---------------------------------------------------------------------------

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
// BamGroupStreamer: reads BAM and yields batches of ReadGroups
// ---------------------------------------------------------------------------

struct BamGroupStreamer<R: Read> {
    reader: noodles_bam::io::Reader<R>,
    batch_size: usize,
    finished: bool,
    paired: bool,
    pending: Option<(String, ReadInfo)>,
    record_buf: noodles_bam::Record,
    total_groups: usize,
    total_alignments: usize,
}

impl<R: Read> BamGroupStreamer<R> {
    fn new(
        reader: noodles_bam::io::Reader<R>,
        batch_size: usize,
    ) -> Self {
        // Read header (already consumed by caller, so we skip this)
        // Actually the caller must NOT have read the header yet from this reader.
        // But in our design, the caller reads the header, creates references, then
        // passes the reader here. So we don't read header again.
        BamGroupStreamer {
            reader,
            batch_size,
            finished: false,
            paired: false,
            pending: None,
            record_buf: noodles_bam::Record::default(),
            total_groups: 0,
            total_alignments: 0,
        }
    }

    /// Read the next batch of ReadGroups. Returns None when BAM is exhausted.
    fn next_batch(&mut self) -> Option<Vec<ReadGroup>> {
        if self.finished {
            return None;
        }

        let mut groups: Vec<ReadGroup> = Vec::with_capacity(self.batch_size);
        let mut current_name: Option<String> = None;
        let mut current_read1s: Vec<ReadInfo> = Vec::new();
        let mut current_read2s: Vec<ReadInfo> = Vec::new();
        let mut current_singles: Vec<ReadInfo> = Vec::new();

        // Seed with pending record from previous batch
        if let Some((name, info)) = self.pending.take() {
            if info.is_paired { self.paired = true; }
            current_name = Some(name);
            if self.paired {
                if info.is_read1 { current_read1s.push(info); }
                else if info.is_read2 { current_read2s.push(info); }
            } else {
                current_singles.push(info);
            }
        }

        loop {
            // Read next record
            let got = match self.reader.read_record(&mut self.record_buf) {
                Ok(0) => false,
                Ok(_) => true,
                Err(_) => continue,
            };

            if !got {
                // EOF: finalize current group
                self.finalize_group(
                    &mut current_name,
                    &mut current_read1s,
                    &mut current_read2s,
                    &mut current_singles,
                    &mut groups,
                );
                self.finished = true;
                break;
            }

            let parsed = parse_record(&self.record_buf);

            match parsed {
                Some((name, info)) => {
                    self.total_alignments += 1;
                    if info.is_paired { self.paired = true; }

                    // Check if this starts a new group
                    let name_changed = match &current_name {
                        Some(cur) => &name != cur,
                        None => false, // first record, start new group
                    };

                    if name_changed {
                        // Finalize current group
                        self.finalize_group(
                            &mut current_name,
                            &mut current_read1s,
                            &mut current_read2s,
                            &mut current_singles,
                            &mut groups,
                        );

                        // Batch full? Save this record as pending and return.
                        if groups.len() >= self.batch_size {
                            self.pending = Some((name, info));
                            break;
                        }
                    }

                    // Add to current group
                    current_name = Some(name);
                    if self.paired {
                        if info.is_read1 { current_read1s.push(info); }
                        else if info.is_read2 { current_read2s.push(info); }
                    } else {
                        current_singles.push(info);
                    }
                }
                None => {
                    // Skipped record (unmapped/duplicate/qc-fail) — continue reading
                    continue;
                }
            }
        }

        self.total_groups += groups.len();
        if groups.is_empty() { None } else { Some(groups) }
    }

    fn finalize_group(
        &self,
        current_name: &mut Option<String>,
        current_read1s: &mut Vec<ReadInfo>,
        current_read2s: &mut Vec<ReadInfo>,
        current_singles: &mut Vec<ReadInfo>,
        groups: &mut Vec<ReadGroup>,
    ) {
        if current_name.is_none() {
            return;
        }
        *current_name = None;

        if self.paired {
            // Move data out
            let r1s = std::mem::take(current_read1s);
            let r2s = std::mem::take(current_read2s);
            groups.push(ReadGroup::Paired { read1s: r1s, read2s: r2s });
        } else {
            let reads = std::mem::take(current_singles);
            groups.push(ReadGroup::Single { reads });
        }
    }
}

// ---------------------------------------------------------------------------
// Aggregation helper (shared between streaming and fallback paths)
// ---------------------------------------------------------------------------

struct AggState {
    gene_counts: HashMap<String, f64>,
    te_counts: Vec<f64>,
    te_multi_counts: Vec<f64>,
    multi_reads: Vec<Vec<usize>>,
    leftover_gene: Vec<(Vec<Vec<String>>, f64)>,
    leftover_te: Vec<(Vec<Vec<usize>>, f64)>,
    empty: i64,
    nonunique: i64,
    uniq_reads: i64,
    avg_read_length: i64,
    tmp_cnt: i64,
}

impl AggState {
    fn new(gene_index: &GeneIndex, te_index: &TEIndex) -> Self {
        let mut gene_counts = HashMap::new();
        for f in gene_index.features() {
            gene_counts.insert(f, 0.0);
        }
        let n_te = te_index.num_instances();
        AggState {
            gene_counts,
            te_counts: vec![0.0; n_te],
            te_multi_counts: vec![0.0; n_te],
            multi_reads: Vec::new(),
            leftover_gene: Vec::new(),
            leftover_te: Vec::new(),
            empty: 0,
            nonunique: 0,
            uniq_reads: 0,
            avg_read_length: 0,
            tmp_cnt: 0,
        }
    }

    fn aggregate(&mut self, annot_results: &[AnnotResult], te_mode: &str, max_length: i64) {
        for res in annot_results {
            if res.is_multi {
                self.nonunique += 1;
                if te_mode == "uniq" {
                    self.empty += 1;
                    continue;
                }
            } else {
                self.uniq_reads += 1;
                if let Some(fl) = res.frag_len {
                    if fl <= max_length && self.tmp_cnt < 10000 {
                        self.avg_read_length += fl;
                        self.tmp_cnt += 1;
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
                    &res.annot_te, &mut self.te_counts, &mut self.te_multi_counts,
                    &mut self.multi_reads, &mut self.leftover_te,
                );
                if no_annot_te {
                    let no_annot_gene = parse_annotations_gene(
                        &res.annot_gene, &mut self.gene_counts, &mut self.leftover_gene,
                    );
                    if no_annot_gene { self.empty += 1; }
                }
            } else {
                let no_annot_gene = parse_annotations_gene(
                    &res.annot_gene, &mut self.gene_counts, &mut self.leftover_gene,
                );
                if no_annot_gene {
                    let no_annot_te = parse_annotations_te(
                        &res.annot_te, &mut self.te_counts, &mut self.te_multi_counts,
                        &mut self.multi_reads, &mut self.leftover_te,
                    );
                    if no_annot_te { self.empty += 1; }
                }
            }
        }
    }
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
    if sort_by_pos {
        // Try samtools sort -n to create a name-sorted temp BAM
        match count_with_external_sort(bam_path, gene_index, te_index, stranded,
                                       te_mode, num_iterations, frag_length, max_length) {
            Ok(result) => return result,
            Err(e) => {
                eprintln!("samtools not available ({}), falling back to in-memory sort", e);
                return count_in_memory(bam_path, gene_index, te_index, stranded,
                                       te_mode, num_iterations, frag_length, max_length);
            }
        }
    }

    // Streaming path (BAM already sorted by name)
    let mut reader = noodles_bam::io::reader::Builder::default()
        .build_from_path(bam_path)
        .unwrap_or_else(|e| {
            eprintln!("Error opening BAM file {}: {}", bam_path, e);
            std::process::exit(1);
        });

    let _header = reader.read_header().expect("Error reading BAM header");
    // We need to rebuild references from the header. noodles consumed it.
    // Re-open to get header separately.
    // Actually we can get references from the header before passing reader to streamer.
    // Let's get references first.
    drop(reader);

    // Re-open: read header, get references, then pass reader to streamer
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

    eprintln!("Reading BAM file (streaming mode)...");
    let batch_size = 50_000;
    let mut streamer = BamGroupStreamer::new(reader, batch_size);
    let mut state = AggState::new(gene_index, te_index);
    let mut batch_num = 0;

    while let Some(groups) = streamer.next_batch() {
        batch_num += 1;
        let annot_results: Vec<AnnotResult> = groups
            .par_iter()
            .map(|group| overlap_annotation_readgroup(group, &references, gene_index, te_index, stranded))
            .collect();

        state.aggregate(&annot_results, te_mode, max_length);
        // groups and annot_results dropped here — memory freed
    }

    eprintln!("Processed {} alignments in {} groups across {} batches.",
              streamer.total_alignments, streamer.total_groups, batch_num);

    finish_count(state, te_index, streamer.paired, num_iterations, frag_length)
}

/// Use samtools sort -n to create a name-sorted temp BAM, then stream it.
fn count_with_external_sort(
    bam_path: &str,
    gene_index: &GeneIndex,
    te_index: &TEIndex,
    stranded: &str,
    te_mode: &str,
    num_iterations: i32,
    frag_length: i64,
    max_length: i64,
) -> Result<CountResult, String> {
    // Check samtools availability
    let check = std::process::Command::new("samtools")
        .arg("--version")
        .output()
        .map_err(|_| "samtools not found in PATH".to_string())?;
    if !check.status.success() {
        return Err("samtools --version failed".to_string());
    }

    let tmp_dir = std::env::temp_dir();
    let tmp_bam = tmp_dir.join(format!("tetranscripts_namesort_{}.bam", std::process::id()));
    let tmp_bam_str = tmp_bam.to_string_lossy().to_string();

    eprintln!("Sorting BAM by read name with samtools...");
    let status = std::process::Command::new("samtools")
        .args(["sort", "-n", "-@", "2", "-o", &tmp_bam_str, bam_path])
        .status()
        .map_err(|e| format!("samtools sort failed: {}", e))?;

    if !status.success() {
        let _ = std::fs::remove_file(&tmp_bam);
        return Err("samtools sort -n exited with error".to_string());
    }

    eprintln!("Streaming sorted BAM...");
    let result = count_streaming_from_path(&tmp_bam_str, gene_index, te_index, stranded,
                                            te_mode, num_iterations, frag_length, max_length);

    let _ = std::fs::remove_file(&tmp_bam);
    result
}

/// Stream from a specific BAM path (used after external sort).
fn count_streaming_from_path(
    bam_path: &str,
    gene_index: &GeneIndex,
    te_index: &TEIndex,
    stranded: &str,
    te_mode: &str,
    num_iterations: i32,
    frag_length: i64,
    max_length: i64,
) -> Result<CountResult, String> {
    let mut reader = noodles_bam::io::reader::Builder::default()
        .build_from_path(bam_path)
        .map_err(|e| format!("Error opening {}: {}", bam_path, e))?;

    let header = reader.read_header().map_err(|e| format!("Error reading header: {}", e))?;
    let references: Vec<String> = header.reference_sequences()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    let batch_size = 50_000;
    let mut streamer = BamGroupStreamer::new(reader, batch_size);
    let mut state = AggState::new(gene_index, te_index);

    while let Some(groups) = streamer.next_batch() {
        let annot_results: Vec<AnnotResult> = groups
            .par_iter()
            .map(|group| overlap_annotation_readgroup(group, &references, gene_index, te_index, stranded))
            .collect();
        state.aggregate(&annot_results, te_mode, max_length);
    }

    Ok(finish_count(state, te_index, streamer.paired, num_iterations, frag_length))
}

/// Fallback: load all into memory, sort, group, process (original approach but with SmallVec).
fn count_in_memory(
    bam_path: &str,
    gene_index: &GeneIndex,
    te_index: &TEIndex,
    stranded: &str,
    te_mode: &str,
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

    let header = reader.read_header().expect("Error reading BAM header");
    let references: Vec<String> = header.reference_sequences()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    eprintln!("Reading BAM file (in-memory fallback)...");
    let mut names: Vec<String> = Vec::new();
    let mut infos: Vec<ReadInfo> = Vec::new();
    let mut record = noodles_bam::Record::default();
    let mut paired = false;

    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {}
            Err(_) => continue,
        }
        if let Some((name, info)) = parse_record(&record) {
            if info.is_paired { paired = true; }
            names.push(name);
            infos.push(info);
        }
    }
    eprintln!("Read {} alignments.", names.len());

    // Sort by name (parallel)
    eprintln!("Sorting by read name...");
    let mut indices: Vec<usize> = (0..names.len()).collect();
    indices.par_sort_by(|&a, &b| names[a].cmp(&names[b]));

    // Group by name and process in batches
    let batch_size = 50_000;
    let mut state = AggState::new(gene_index, te_index);
    let mut batch_groups: Vec<ReadGroup> = Vec::with_capacity(batch_size);
    let mut i = 0;
    let mut batch_num = 0;

    while i < indices.len() {
        let start = i;
        let name = &names[indices[start]];

        while i < indices.len() && &names[indices[i]] == name {
            i += 1;
        }

        // Build group
        if paired {
            let mut read1s = Vec::new();
            let mut read2s = Vec::new();
            for j in start..i {
                let info = infos[indices[j]].clone();
                if info.is_read1 { read1s.push(info.clone()); }
                if info.is_read2 { read2s.push(info); }
            }
            batch_groups.push(ReadGroup::Paired { read1s, read2s });
        } else {
            let reads: Vec<ReadInfo> = (start..i).map(|j| infos[indices[j]].clone()).collect();
            batch_groups.push(ReadGroup::Single { reads });
        }

        if batch_groups.len() >= batch_size {
            batch_num += 1;
            let annot_results: Vec<AnnotResult> = batch_groups
                .par_iter()
                .map(|group| overlap_annotation_readgroup(group, &references, gene_index, te_index, stranded))
                .collect();
            state.aggregate(&annot_results, te_mode, max_length);
            batch_groups.clear();
        }
    }

    // Process remaining
    if !batch_groups.is_empty() {
        batch_num += 1;
        let annot_results: Vec<AnnotResult> = batch_groups
            .par_iter()
            .map(|group| overlap_annotation_readgroup(group, &references, gene_index, te_index, stranded))
            .collect();
        state.aggregate(&annot_results, te_mode, max_length);
    }

    eprintln!("Processed {} batches.", batch_num);
    finish_count(state, te_index, paired, num_iterations, frag_length)
}

// ---------------------------------------------------------------------------
// Finalize: resolve ambiguity, EM, return CountResult
// ---------------------------------------------------------------------------

fn finish_count(
    mut state: AggState,
    te_index: &TEIndex,
    paired: bool,
    num_iterations: i32,
    frag_length: i64,
) -> CountResult {
    // Resolve leftover ambiguities
    if !state.leftover_gene.is_empty() {
        resolve_annotation_ambiguity_gene(&mut state.gene_counts, &state.leftover_gene);
    }
    if !state.leftover_te.is_empty() {
        resolve_annotation_ambiguity_te(&mut state.te_counts, &state.leftover_te);
    }

    eprintln!("uniq te counts = {}", state.te_counts.iter().sum::<f64>() as i64);

    let estimated_read_length = if !paired && frag_length > 0 {
        frag_length
    } else if state.avg_read_length > 0 && state.tmp_cnt > 0 {
        state.avg_read_length / state.tmp_cnt
    } else {
        100
    };

    let new_te_multi_counts = if num_iterations > 0 && !state.multi_reads.is_empty() {
        eprintln!("Starting iterative optimization...");
        em_estimate(te_index, &state.multi_reads, &state.te_counts, &state.te_multi_counts,
                    num_iterations, estimated_read_length)
    } else {
        state.te_multi_counts
    };

    for j in 0..state.te_counts.len() {
        state.te_counts[j] += new_te_multi_counts[j];
    }

    let st: f64 = state.te_counts.iter().sum();
    let sg: f64 = state.gene_counts.values().sum();

    eprintln!("TE counts total {}", st);
    eprintln!("Gene counts total {}", sg);
    eprintln!("Total annotated = {}", (st + sg) as i64);
    eprintln!("Total non-unique = {}", state.nonunique);
    eprintln!("Total unannotated = {}", state.empty);

    CountResult {
        gene_counts: state.gene_counts,
        te_instance_counts: state.te_counts,
        total_annotated: (st + sg) as i64,
        total_nonunique: state.nonunique,
        total_unannotated: state.empty,
    }
}
