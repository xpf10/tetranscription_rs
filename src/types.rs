use pyo3::prelude::*;

/// Strand direction for genomic intervals
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Strand {
    pub fn from_char(s: &str) -> Self {
        match s {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => Strand::Unknown,
        }
    }

    pub fn to_direction(&self) -> i32 {
        match self {
            Strand::Plus => 1,
            Strand::Minus => -1,
            Strand::Unknown => 0,
        }
    }
}

// Implement IntoPyObject for Strand so PyO3 can convert it
impl<'py> IntoPyObject<'py> for Strand {
    type Target = pyo3::types::PyString;
    type Output = pyo3::Bound<'py, Self::Target>;
    type Error = std::convert::Infallible;

    fn into_pyobject(self, py: pyo3::Python<'py>) -> Result<Self::Output, Self::Error> {
        Ok(pyo3::types::PyString::new(py, match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
            Strand::Unknown => ".",
        }))
    }
}

/// A genomic interval representing an exon from a read alignment
#[derive(Debug, Clone)]
pub struct ExonInterval {
    pub chrom: String,
    pub start: i64,  // 1-based
    pub end: i64,
    pub strand: Strand,
}

/// CIGAR operation code and length
#[derive(Debug, Clone, Copy)]
pub struct CigarElement {
    pub code: u32,
    pub len: i64,
}

/// Constants
pub const OPT_TOL: f64 = 0.0001;
