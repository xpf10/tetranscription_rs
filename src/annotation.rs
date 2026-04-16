/// Annotation parsing and ambiguity resolution.
/// Ported from bin/TEcount: parse_annotations_gene, parse_annotations_TE, resolve_annotation_ambiguity

use std::collections::HashMap;

/// Parse gene annotations and update gene_counts.
/// Returns true if there was no annotation (unannotated read).
pub fn parse_annotations_gene(
    annot_gene: &[Vec<String>],
    gene_counts: &mut HashMap<String, f64>,
    leftover_gene: &mut Vec<(Vec<Vec<String>>, f64)>,
) -> bool {
    if annot_gene.len() > 1 {
        leftover_gene.push((annot_gene.to_vec(), 1.0));
    } else if annot_gene.len() == 1 {
        let genes = &annot_gene[0];
        if genes.len() == 1 {
            *gene_counts.entry(genes[0].clone()).or_insert(0.0) += 1.0;
        } else {
            let w = 1.0 / genes.len() as f64;
            for g in genes {
                *gene_counts.entry(g.clone()).or_insert(0.0) += w;
            }
        }
    } else {
        return true; // no annotation
    }
    false
}

/// Parse TE annotations and update counts.
/// Returns true if there was no annotation (unannotated read).
pub fn parse_annotations_te(
    annot_te: &[Vec<usize>],
    te_counts: &mut Vec<f64>,
    te_multi_counts: &mut Vec<f64>,
    multi_reads: &mut Vec<Vec<usize>>,
    leftover_te: &mut Vec<(Vec<Vec<usize>>, f64)>,
) -> bool {
    if annot_te.is_empty() {
        return true;
    }

    if annot_te.len() == 1 && annot_te[0].len() == 1 {
        let idx = annot_te[0][0];
        if idx < te_counts.len() {
            te_counts[idx] += 1.0;
        }
    }

    if annot_te.len() == 1 && annot_te[0].len() > 1 {
        leftover_te.push((annot_te.to_vec(), 1.0));
    }

    if annot_te.len() > 1 {
        let mut multi_algn = Vec::new();
        for annot in annot_te {
            for &te in annot {
                let w = 1.0 / (annot_te.len() * annot.len()) as f64;
                if te < te_multi_counts.len() {
                    te_multi_counts[te] += w;
                }
                multi_algn.push(te);
            }
        }
        multi_reads.push(multi_algn);
    }

    false
}

/// Resolve annotation ambiguity for TE counts (weighted assignment).
pub fn resolve_annotation_ambiguity_te(
    counts: &mut Vec<f64>,
    leftovers: &[(Vec<Vec<usize>>, f64)],
) {
    for (annlist, w) in leftovers {
        let mut readslist: HashMap<usize, f64> = HashMap::new();
        let mut total = 0.0;
        let size = annlist.len();
        let ww = if size > 1 { w / size as f64 } else { *w };

        for ann in annlist {
            for &a in ann {
                let entry = counts.get(a).copied().unwrap_or(0.0);
                let e = readslist.entry(a).or_insert(0.0);
                if *e == 0.0 {
                    *e = entry;
                    total += entry;
                }
            }
        }

        if total > 0.0 {
            for (&a, &v) in &readslist {
                let add = ww * v / total;
                if a < counts.len() {
                    counts[a] += add;
                }
            }
        } else {
            let n = readslist.len() as f64;
            for &a in readslist.keys() {
                if a < counts.len() {
                    counts[a] += ww / n;
                }
            }
        }
    }
}

/// Resolve annotation ambiguity for gene counts (weighted assignment).
pub fn resolve_annotation_ambiguity_gene(
    counts: &mut HashMap<String, f64>,
    leftovers: &[(Vec<Vec<String>>, f64)],
) {
    for (annlist, w) in leftovers {
        let mut readslist: HashMap<String, f64> = HashMap::new();
        let mut total = 0.0;
        let size = annlist.len();
        let ww = if size > 1 { w / size as f64 } else { *w };

        for ann in annlist {
            for a in ann {
                let entry = counts.get(a).copied().unwrap_or(0.0);
                let e = readslist.entry(a.clone()).or_insert(0.0);
                if *e == 0.0 {
                    *e = entry;
                    total += entry;
                }
            }
        }

        if total > 0.0 {
            for (a, v) in &readslist {
                let add = ww * v / total;
                *counts.entry(a.clone()).or_insert(0.0) += add;
            }
        } else {
            let n = readslist.len() as f64;
            for a in readslist.keys() {
                *counts.entry(a.clone()).or_insert(0.0) += ww / n;
            }
        }
    }
}
