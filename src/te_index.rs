/// TE annotation index using BTreeMap with bin-based bucketing.
/// Ported from TEToolkit/TEindex.py (TEfeatures, BinaryTree, Node)
/// Replaced custom AVL tree with std::collections::BTreeMap for balanced O(log n) operations.

use std::collections::{BTreeMap, HashMap, HashSet};

use pyo3::prelude::*;

use crate::gtf_parser::parse_te_gtf;
use crate::types::{ExonInterval, Strand, TEINDEX_BINSIZE};

// ---------------------------------------------------------------------------
// Bin-entry: stores (name_idx, end) pairs keyed by start position
// ---------------------------------------------------------------------------

struct BinEntries {
    entries: HashMap<i64, Vec<(usize, i64)>>,
}

impl BinEntries {
    fn new(start: i64, end: i64, name_idx: usize) -> Self {
        let mut entries = HashMap::new();
        entries.insert(start, vec![(name_idx, end)]);
        BinEntries { entries }
    }

    fn add(&mut self, start: i64, end: i64, name_idx: usize) {
        self.entries.entry(start).or_default().push((name_idx, end));
    }

    fn overlaps(&self, start: i64, end: i64) -> Vec<usize> {
        let mut result = Vec::new();
        for (&s, pairs) in &self.entries {
            if s > end {
                continue;
            }
            for &(idx, e) in pairs {
                if start <= e && end >= s {
                    result.push(idx);
                }
            }
        }
        result
    }
}

// ---------------------------------------------------------------------------
// TE tree: BTreeMap<bin_start, BinEntries> per chromosome
// ---------------------------------------------------------------------------

struct TeTree {
    bins: BTreeMap<i64, BinEntries>,
}

impl TeTree {
    fn new() -> Self {
        TeTree { bins: BTreeMap::new() }
    }

    fn insert(&mut self, start: i64, end: i64, name_idx: usize) {
        let bin_start = bin_start_id(start);
        self.bins
            .entry(bin_start)
            .and_modify(|e| e.add(start, end, name_idx))
            .or_insert_with(|| BinEntries::new(start, end, name_idx));
    }

    /// Find all TE instances overlapping the genomic interval [query_start, query_end].
    /// Uses bin range for efficient lookup, then checks actual genomic coordinates.
    fn find_overlapping(&self, query_start: i64, query_end: i64) -> Vec<usize> {
        let start_bin = bin_start_id(query_start);
        let end_bin = bin_end_id(query_end);
        let mut result = Vec::new();
        for (_, entries) in self.bins.range(start_bin..=end_bin) {
            result.extend(entries.overlaps(query_start, query_end));
        }
        result
    }
}

fn bin_start_id(pos: i64) -> i64 {
    let mut bid = pos / TEINDEX_BINSIZE;
    if pos == bid * TEINDEX_BINSIZE {
        bid -= 1;
    }
    bid
}

#[allow(dead_code)]
fn bin_end_id(pos: i64) -> i64 {
    pos / TEINDEX_BINSIZE
}

// ---------------------------------------------------------------------------
// TEIndex
// ---------------------------------------------------------------------------

#[pyclass]
pub struct TEIndex {
    trees: HashMap<String, TeTree>,
    lengths: Vec<i64>,
    name_id_map: Vec<String>,
    elements: Vec<String>,
}

#[pymethods]
impl TEIndex {
    #[new]
    pub fn from_gtf(tefile: &str) -> Self {
        let records = parse_te_gtf(tefile);
        let mut trees: HashMap<String, TeTree> = HashMap::new();
        let mut lengths = Vec::new();
        let mut name_id_map = Vec::new();
        let mut elements: Vec<String> = Vec::new();
        let mut seen_elements: HashSet<String> = HashSet::new();

        for rec in &records {
            let tlen = rec.end - rec.start + 1;

            let transcript_id = rec.transcript_id.as_deref().unwrap_or("");
            let family_id = rec.family_id.as_deref().unwrap_or("");
            let class_id = rec.class_id.as_deref().unwrap_or("");
            let strand = &rec.strand;

            let full_name = format!("{}:{}:{}:{}:{}", transcript_id, rec.gene_id, family_id, class_id, strand);
            let ele_name = format!("{}:{}:{}", rec.gene_id, family_id, class_id);

            if !seen_elements.contains(&ele_name) {
                seen_elements.insert(ele_name.clone());
                elements.push(ele_name);
            }

            let name_idx = name_id_map.len();
            lengths.push(tlen);
            name_id_map.push(full_name);

            let tree = trees.entry(rec.chrom.clone()).or_insert_with(TeTree::new);

            let mut bin_start = bin_start_id(rec.start);
            let bin_end = bin_end_id(rec.end);

            while bin_start <= bin_end {
                let end_pos = std::cmp::min(rec.end, (bin_start + 1) * TEINDEX_BINSIZE);
                let start_pos = std::cmp::max(rec.start, bin_start * TEINDEX_BINSIZE + 1);
                tree.insert(start_pos, end_pos, name_idx);
                bin_start += 1;
            }
        }

        TEIndex { trees, lengths, name_id_map, elements }
    }

    pub fn num_instances(&self) -> usize {
        self.name_id_map.len()
    }

    pub fn get_length(&self, idx: usize) -> i64 {
        if idx < self.lengths.len() { self.lengths[idx] } else { -1 }
    }

    pub fn get_ele_name(&self, idx: usize) -> Option<String> {
        if idx >= self.name_id_map.len() { return None; }
        let full_name = &self.name_id_map[idx];
        let parts: Vec<&str> = full_name.split(':').collect();
        if parts.len() >= 4 {
            Some(format!("{}:{}:{}", parts[1], parts[2], parts[3]))
        } else {
            None
        }
    }

    pub fn get_strand(&self, idx: usize) -> Strand {
        if idx >= self.name_id_map.len() { return Strand::Unknown; }
        let full_name = &self.name_id_map[idx];
        if let Some(c) = full_name.chars().last() {
            match c {
                '+' => Strand::Plus,
                '-' => Strand::Minus,
                _ => Strand::Unknown,
            }
        } else {
            Strand::Unknown
        }
    }

    pub fn group_by_element(&self, instance_counts: Vec<f64>) -> HashMap<String, f64> {
        let mut ele_counts: HashMap<String, f64> = HashMap::new();
        for ele in &self.elements {
            ele_counts.insert(ele.clone(), 0.0);
        }
        for (i, &cnt) in instance_counts.iter().enumerate() {
            if let Some(ele_name) = self.get_ele_name(i) {
                *ele_counts.entry(ele_name).or_insert(0.0) += cnt;
            }
        }
        ele_counts
    }

    pub fn get_elements(&self) -> Vec<String> {
        self.elements.clone()
    }
}

impl TEIndex {
    pub fn te_annotation(&self, intervals: &[ExonInterval]) -> Vec<usize> {
        let mut tes = Vec::new();
        let mut seen: HashSet<usize> = HashSet::new();

        for itv in intervals {
            if let Some(name_idx_list) = self.find_overlapping_te(&itv.chrom, itv.start, itv.end) {
                for t in name_idx_list {
                    let te_strand = self.get_strand(t);
                    let strand_match = match itv.strand {
                        Strand::Unknown => true,
                        _ => te_strand == itv.strand,
                    };
                    if strand_match && !seen.contains(&t) {
                        seen.insert(t);
                        tes.push(t);
                    }
                }
            }
        }
        tes
    }

    fn find_overlapping_te(&self, chrom: &str, start: i64, end: i64) -> Option<Vec<usize>> {
        let tree = self.trees.get(chrom)?;
        Some(tree.find_overlapping(start, end))
    }
}
