/// GTF file parser for gene and TE annotations.
/// Ported from TEToolkit/GeneFeatures.py (GFF_Reader) and TEToolkit/TEindex.py (TEfeatures.build)

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// A parsed GTF record
#[derive(Debug, Clone)]
pub struct GtfRecord {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
    pub feature: String,
    pub gene_id: String,
    pub transcript_id: Option<String>,
    pub family_id: Option<String>,
    pub class_id: Option<String>,
}

/// Parse all attributes from a GTF attribute string
fn parse_all_attrs(attr_str: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    for pair in attr_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        if let Some(pos) = pair.find(' ') {
            let k = pair[..pos].trim().to_string();
            let mut v = pair[pos + 1..].trim().to_string();
            v = v.replace('"', "");
            attrs.insert(k, v);
        }
    }
    attrs
}

/// Parse a gene GTF file. Returns (records, gene_id_list).
pub fn parse_gene_gtf(
    path: &str,
    feature_type: &str,
    id_attribute: &str,
) -> (Vec<GtfRecord>, Vec<String>) {
    let mut records = Vec::new();
    let mut gene_ids: Vec<String> = Vec::new();
    let mut seen_genes = std::collections::HashSet::new();

    let file = File::open(path).unwrap_or_else(|e| {
        eprintln!("Error opening gene GTF file {}: {}", path, e);
        std::process::exit(1);
    });
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let feature = fields[2];
        if feature != feature_type {
            continue;
        }

        let attrs = parse_all_attrs(fields[8]);
        let gene_id = match attrs.get(id_attribute) {
            Some(id) => id.clone(),
            None => continue,
        };

        let record = GtfRecord {
            chrom: fields[0].to_string(),
            start: fields[3].parse::<i64>().unwrap_or(0),
            end: fields[4].parse::<i64>().unwrap_or(0),
            strand: fields[6].to_string(),
            feature: feature.to_string(),
            gene_id: gene_id.clone(),
            transcript_id: None,
            family_id: None,
            class_id: None,
        };

        if !seen_genes.contains(&gene_id) {
            seen_genes.insert(gene_id.clone());
            gene_ids.push(gene_id.clone());
        }

        records.push(record);
    }

    (records, gene_ids)
}

/// Parse a TE GTF file. Records are sorted by chromosome and start position.
pub fn parse_te_gtf(path: &str) -> Vec<GtfRecord> {
    let file = File::open(path).unwrap_or_else(|e| {
        eprintln!("Error opening TE GTF file {}: {}", path, e);
        std::process::exit(1);
    });
    let reader = BufReader::new(file);
    let mut raw_lines: Vec<String> = reader.lines().filter_map(|l| l.ok()).collect();

    raw_lines.sort_by(|a, b| {
        let a_fields: Vec<&str> = a.split('\t').collect();
        let b_fields: Vec<&str> = b.split('\t').collect();
        let a_chrom = a_fields.first().unwrap_or(&"");
        let b_chrom = b_fields.first().unwrap_or(&"");
        match a_chrom.cmp(b_chrom) {
            std::cmp::Ordering::Equal => {
                let a_start: i64 = a_fields.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);
                let b_start: i64 = b_fields.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);
                a_start.cmp(&b_start)
            }
            other => other,
        }
    });

    let mut records = Vec::new();
    for line in &raw_lines {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let attrs = parse_all_attrs(fields[8]);
        let gene_id = match attrs.get("gene_id") {
            Some(id) => id.clone(),
            None => continue,
        };

        let record = GtfRecord {
            chrom: fields[0].to_string(),
            start: fields[3].parse::<i64>().unwrap_or(0),
            end: fields[4].parse::<i64>().unwrap_or(0),
            strand: fields[6].to_string(),
            feature: fields[2].to_string(),
            gene_id,
            transcript_id: attrs.get("transcript_id").cloned(),
            family_id: attrs.get("family_id").cloned(),
            class_id: attrs.get("class_id").cloned(),
        };

        records.push(record);
    }

    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_all_attrs() {
        let attr = r#"gene_id "ENSG00001"; transcript_id "ENST00001";"#;
        let attrs = parse_all_attrs(attr);
        assert_eq!(attrs.get("gene_id"), Some(&"ENSG00001".to_string()));
        assert_eq!(attrs.get("transcript_id"), Some(&"ENST00001".to_string()));
    }
}
