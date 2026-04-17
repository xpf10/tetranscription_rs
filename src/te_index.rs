/// TE annotation index using flat sorted arrays with end-max augmentation.
/// Replaces the previous BTreeMap<HashMap<Vec>> design to reduce memory
/// from ~2.6 GB to ~120 MB for 3.7M TE instances.

use std::collections::{HashMap, HashSet};

use pyo3::prelude::*;

use crate::gtf_parser::parse_te_gtf;
use crate::types::{ExonInterval, Strand};

// ---------------------------------------------------------------------------
// TeTree: flat sorted array with end-max augmentation for O(log n + k) queries
// ---------------------------------------------------------------------------

struct TeTree {
    /// Intervals sorted by start position: (start, end, name_idx)
    intervals: Vec<(i64, i64, usize)>,
    /// max_ends[i] = max(intervals[i..].1), used for pruning backward scan
    max_ends: Vec<i64>,
}

impl TeTree {
    fn new(mut entries: Vec<(i64, i64, usize)>) -> Self {
        // Sort by start position
        entries.sort_by_key(|(s, _, _)| *s);

        // Build max_ends: scan from right to left
        let n = entries.len();
        let mut max_ends = vec![0i64; n];
        if n > 0 {
            let mut max_end = entries[n - 1].1;
            for i in (0..n).rev() {
                if entries[i].1 > max_end {
                    max_end = entries[i].1;
                }
                max_ends[i] = max_end;
            }
        }

        TeTree {
            intervals: entries,
            max_ends,
        }
    }

    /// Find all TE instances overlapping [query_start, query_end].
    /// Uses binary search + end-max pruning for O(log n + k) queries.
    fn find_overlapping(&self, query_start: i64, query_end: i64) -> Vec<usize> {
        if self.intervals.is_empty() {
            return Vec::new();
        }

        // Binary search: find the last interval with start <= query_end
        let end_idx = match self.intervals.binary_search_by_key(&query_end, |(s, _, _)| *s) {
            Ok(idx) => idx + 1,
            Err(idx) => idx,
        };

        if end_idx == 0 {
            return Vec::new();
        }

        let mut result = Vec::new();
        let mut i = end_idx;
        while i > 0 {
            i -= 1;
            // Pruning: if no interval from 0..=i can reach into [query_start, query_end]
            if self.max_ends[i] < query_start {
                break;
            }
            let (start, end, name_idx) = self.intervals[i];
            if end >= query_start && start <= query_end {
                result.push(name_idx);
            }
        }
        result
    }
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
        let mut tree_entries: HashMap<String, Vec<(i64, i64, usize)>> = HashMap::new();
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

            // Store the TE interval directly — no bin splitting needed
            tree_entries
                .entry(rec.chrom.clone())
                .or_default()
                .push((rec.start, rec.end, name_idx));
        }

        // Build TeTrees from sorted entries
        let trees: HashMap<String, TeTree> = tree_entries
            .into_iter()
            .map(|(chrom, entries)| (chrom, TeTree::new(entries)))
            .collect();

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
