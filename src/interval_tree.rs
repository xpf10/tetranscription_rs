/// Center-based interval tree for gene annotation overlap queries.
/// Ported from TEToolkit/IntervalTree.py

use std::cmp::Ordering;

/// A genomic interval associated with a gene ID
#[derive(Debug, Clone)]
pub struct Interval {
    pub gene_id: String,
    pub start: i64,
    pub stop: i64,
}

impl Interval {
    pub fn new(gene_id: String, start: i64, stop: i64) -> Self {
        Interval { gene_id, start, stop }
    }
}

impl PartialEq for Interval {
    fn eq(&self, other: &Self) -> bool {
        self.start == other.start && self.stop == other.stop && self.gene_id == other.gene_id
    }
}

impl Eq for Interval {}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> Ordering {
        self.start.cmp(&other.start).then(self.stop.cmp(&other.stop))
    }
}

const DEFAULT_DEPTH: usize = 16;
const DEFAULT_MIN_BUCKET: usize = 16;
const DEFAULT_MAX_BUCKET: usize = 512;

/// Center-based interval tree. Ported from the original Python implementation.
pub struct IntervalTree {
    intervals: Vec<Interval>,
    left: Option<Box<IntervalTree>>,
    right: Option<Box<IntervalTree>>,
    center: f64,
}

impl IntervalTree {
    /// Build an interval tree from a list of intervals.
    pub fn new(intervals: Vec<Interval>) -> Self {
        Self::build(intervals, DEFAULT_DEPTH, DEFAULT_MIN_BUCKET, None, DEFAULT_MAX_BUCKET)
    }

    fn build(
        mut intervals: Vec<Interval>,
        depth: usize,
        min_bucket: usize,
        extent: Option<(i64, i64)>,
        max_bucket: usize,
    ) -> Self {
        if depth == 0 || (intervals.len() < min_bucket && intervals.len() < max_bucket) {
            return IntervalTree {
                intervals,
                left: None,
                right: None,
                center: 0.0,
            };
        }

        if extent.is_none() {
            intervals.sort();
        }

        let (left_bound, right_bound) = extent.unwrap_or_else(|| {
            let l = intervals[0].start;
            let r = intervals.iter().map(|i| i.stop).max().unwrap_or(0);
            (l, r)
        });

        let center = (left_bound as f64 + right_bound as f64) / 2.0;

        let mut center_intervals = Vec::new();
        let mut lefts = Vec::new();
        let mut rights = Vec::new();

        for interval in intervals {
            if (interval.stop as f64) < center {
                lefts.push(interval);
            } else if (interval.start as f64) > center {
                rights.push(interval);
            } else {
                center_intervals.push(interval);
            }
        }

        let left_tree = if !lefts.is_empty() {
            let left_start = if let Some(ext) = extent {
                ext.0
            } else {
                lefts[0].start
            };
            Some(Box::new(Self::build(
                lefts,
                depth - 1,
                min_bucket,
                Some((left_start, center as i64)),
                max_bucket,
            )))
        } else {
            None
        };

        let right_tree = if !rights.is_empty() {
            let right_end = if let Some(ext) = extent {
                ext.1
            } else {
                rights.iter().map(|i| i.stop).max().unwrap_or(0)
            };
            Some(Box::new(Self::build(
                rights,
                depth - 1,
                min_bucket,
                Some((center as i64, right_end)),
                max_bucket,
            )))
        } else {
            None
        };

        IntervalTree {
            intervals: center_intervals,
            left: left_tree,
            right: right_tree,
            center,
        }
    }

    /// Find all intervals overlapping [start, stop]
    pub fn find(&self, start: i64, stop: i64) -> Vec<&Interval> {
        let mut result = Vec::new();

        if !self.intervals.is_empty() {
            // Check if stop could overlap with our intervals
            if let Some(first) = self.intervals.first() {
                if stop >= first.start {
                    for iv in &self.intervals {
                        if iv.stop >= start && iv.start <= stop {
                            result.push(iv);
                        }
                    }
                }
            }
        }

        if let Some(ref left) = self.left {
            if (start as f64) <= self.center {
                result.extend(left.find(start, stop));
            }
        }

        if let Some(ref right) = self.right {
            if (stop as f64) >= self.center {
                result.extend(right.find(start, stop));
            }
        }

        result
    }

    /// Find gene IDs of all intervals overlapping [start, stop]
    pub fn find_gene(&self, start: i64, stop: i64) -> Vec<String> {
        self.find(start, stop)
            .into_iter()
            .map(|iv| iv.gene_id.clone())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_find() {
        let intervals = vec![
            Interval::new("geneA".into(), 100, 200),
            Interval::new("geneB".into(), 150, 300),
            Interval::new("geneC".into(), 400, 500),
        ];
        let tree = IntervalTree::new(intervals);
        let result = tree.find_gene(120, 180);
        assert!(result.contains(&"geneA".to_string()));
        assert!(result.contains(&"geneB".to_string()));
        assert!(!result.contains(&"geneC".to_string()));
    }

    #[test]
    fn test_no_overlap() {
        let intervals = vec![
            Interval::new("geneA".into(), 100, 200),
            Interval::new("geneB".into(), 300, 400),
        ];
        let tree = IntervalTree::new(intervals);
        let result = tree.find_gene(250, 260);
        assert!(result.is_empty());
    }

    #[test]
    fn test_find_gene_empty() {
        let tree = IntervalTree::new(vec![]);
        let result = tree.find_gene(0, 100);
        assert!(result.is_empty());
    }
}
