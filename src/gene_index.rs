/// Gene annotation index using interval trees.
/// Ported from TEToolkit/GeneFeatures.py (GeneFeatures class)

use std::collections::HashMap;

use pyo3::prelude::*;

use crate::gtf_parser::parse_gene_gtf;
use crate::interval_tree::{Interval, IntervalTree};
use crate::types::ExonInterval;

/// Gene annotation index built from a GTF file.
/// Internally uses three sets of interval trees (plus/minus/no_strand) per chromosome.
#[pyclass]
pub struct GeneIndex {
    features: Vec<String>,
    plus: HashMap<String, IntervalTree>,
    minus: HashMap<String, IntervalTree>,
    no_strand: HashMap<String, IntervalTree>,
}

#[pymethods]
impl GeneIndex {
    /// Build gene index from a GTF file.
    /// - `gtf_path`: path to gene annotation GTF
    /// - `stranded`: "no", "forward", or "reverse"
    /// - `feature_type`: feature to extract (e.g. "exon")
    /// - `id_attribute`: attribute key for gene ID (e.g. "gene_id")
    #[new]
    pub fn from_gtf(
        gtf_path: &str,
        _stranded: &str,
        feature_type: &str,
        id_attribute: &str,
    ) -> Self {
        let (records, features) = parse_gene_gtf(gtf_path, feature_type, id_attribute);

        // Group intervals by chromosome and strand
        let mut temp_plus: HashMap<String, Vec<Interval>> = HashMap::new();
        let mut temp_minus: HashMap<String, Vec<Interval>> = HashMap::new();
        let mut temp_nostrand: HashMap<String, Vec<Interval>> = HashMap::new();

        for rec in &records {
            let iv = Interval::new(rec.gene_id.clone(), rec.start, rec.end);
            match rec.strand.as_str() {
                "+" => {
                    temp_plus.entry(rec.chrom.clone()).or_default().push(iv);
                }
                "-" => {
                    temp_minus.entry(rec.chrom.clone()).or_default().push(iv);
                }
                _ => {
                    temp_nostrand.entry(rec.chrom.clone()).or_default().push(iv);
                }
            }
        }

        // Build interval trees
        let plus = build_tree_map(temp_plus);
        let minus = build_tree_map(temp_minus);
        let no_strand = build_tree_map(temp_nostrand);

        GeneIndex {
            features,
            plus,
            minus,
            no_strand,
        }
    }

    /// Get list of all gene IDs in this index
    #[getter]
    pub fn features(&self) -> Vec<String> {
        self.features.clone()
    }
}

impl GeneIndex {
    /// Find genes overlapping a list of exon intervals.
    /// Returns a list of gene IDs.
    pub fn gene_annotation(&self, itv_list: &[ExonInterval]) -> Vec<String> {
        let mut genes = Vec::new();
        for itv in itv_list {
            let chrom = &itv.chrom;
            let start = itv.start;
            let end = itv.end;

            match itv.strand {
                crate::types::Strand::Plus => {
                    if let Some(tree) = self.plus.get(chrom) {
                        genes.extend(tree.find_gene(start, end));
                    }
                }
                crate::types::Strand::Minus => {
                    if let Some(tree) = self.minus.get(chrom) {
                        genes.extend(tree.find_gene(start, end));
                    }
                }
                crate::types::Strand::Unknown => {
                    // "." strand: check both plus, minus, and no_strand
                    if let Some(tree) = self.minus.get(chrom) {
                        genes.extend(tree.find_gene(start, end));
                    }
                    if let Some(tree) = self.plus.get(chrom) {
                        genes.extend(tree.find_gene(start, end));
                    }
                    if let Some(tree) = self.no_strand.get(chrom) {
                        genes.extend(tree.find_gene(start, end));
                    }
                }
            }
        }
        genes
    }
}

fn build_tree_map(groups: HashMap<String, Vec<Interval>>) -> HashMap<String, IntervalTree> {
    groups
        .into_iter()
        .map(|(chrom, intervals)| (chrom, IntervalTree::new(intervals)))
        .collect()
}
