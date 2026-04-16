use pyo3::prelude::*;

mod annotation;
mod bam_reader;
mod em_algorithm;
mod gene_index;
mod gtf_parser;
mod interval_tree;
mod te_index;
mod types;

use gene_index::GeneIndex;
use te_index::TEIndex;
use bam_reader::count_transcript_abundance_py;

#[pymodule]
fn _core(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GeneIndex>()?;
    m.add_class::<TEIndex>()?;
    m.add_function(wrap_pyfunction!(count_transcript_abundance_py, m)?)?;
    Ok(())
}
