/// EM/SQUAREM iterative optimization for multi-read assignment.
/// Ported from TEToolkit/EMAlgorithm.py

use crate::te_index::TEIndex;
use crate::types::OPT_TOL;

/// Normalize a vector of means to sum to 1.0
fn normalize_means(means: &[f64]) -> Vec<f64> {
    let total: f64 = means.iter().sum();
    if total > 0.0 {
        means.iter().map(|x| x / total).collect()
    } else {
        vec![0.0; means.len()]
    }
}

/// Compute abundances by distributing multi-reads proportionally
fn compute_abundances(means: &[f64], multi_reads: &[Vec<usize>]) -> Vec<f64> {
    let size = means.len();
    let mut multi_counts = vec![0.0; size];

    for te_transcripts in multi_reads {
        let total_mass: f64 = te_transcripts.iter().map(|&tid| means[tid]).sum();
        let norm = if total_mass > 0.0 { 1.0 / total_mass } else { 0.0 };

        for &tid in te_transcripts {
            if tid < size {
                multi_counts[tid] += means[tid] * norm;
            }
        }
    }

    multi_counts
}

/// Single EM update step
fn em_update(
    means_in: &[f64],
    te_index: &TEIndex,
    uniq_counts: &[f64],
    multi_reads: &[Vec<usize>],
    estimated_read_length: i64,
) -> Vec<f64> {
    let multi_counts = compute_abundances(means_in, multi_reads);
    let mut means_out = vec![0.0; means_in.len()];

    for tid in 0..means_in.len() {
        let tlen = te_index.get_length(tid);
        if tlen < 0 {
            eprintln!("Error in optimization: TE {} does not exist!", tid);
            continue;
        }
        let effective_length = tlen - estimated_read_length + 1;
        if effective_length > 0 {
            means_out[tid] = (uniq_counts[tid] + multi_counts[tid]) / effective_length as f64;
        } else {
            means_out[tid] = 0.0;
        }
    }

    normalize_means(&means_out)
}

/// Dot product of two vectors
fn dot_product(u: &[f64], v: &[f64]) -> f64 {
    u.iter().zip(v.iter()).map(|(a, b)| a * b).sum()
}

/// SQUAREM-accelerated EM estimation for multi-read assignment.
/// Returns the estimated multi-read counts.
pub fn em_estimate(
    te_index: &TEIndex,
    multi_reads: &[Vec<usize>],
    uniq_counts: &[f64],
    multi_counts: &[f64],
    num_iterations: i32,
    estimated_read_length: i64,
) -> Vec<f64> {
    let t_size = uniq_counts.len();
    if t_size == 0 || multi_reads.is_empty() {
        return vec![0.0; t_size];
    }

    // Initialize means0: density per base
    let mut means0: Vec<f64> = Vec::with_capacity(t_size);
    for tid in 0..t_size {
        let tlen = te_index.get_length(tid);
        let effective_length = tlen - estimated_read_length + 1;
        if effective_length > 0 {
            means0.push((uniq_counts[tid] + multi_counts[tid]) / effective_length as f64);
        } else {
            means0.push(0.0);
        }
    }
    means0 = normalize_means(&means0);

    let m_step: f64 = 4.0;
    let mut max_step: f64 = 1.0;
    let mut min_step: f64 = 1.0;
    let max_step0: f64 = 1.0;

    let mut cur_iter = 0;
    while cur_iter < num_iterations {
        cur_iter += 1;

        // First EM step
        let means1 = em_update(&means0, te_index, uniq_counts, multi_reads, estimated_read_length);

        // Second EM step
        let means2 = em_update(&means1, te_index, uniq_counts, multi_reads, estimated_read_length);

        // Compute r and v vectors
        let r: Vec<f64> = means1.iter().zip(&means0).map(|(a, b)| a - b).collect();
        let v: Vec<f64> = means2
            .iter()
            .zip(&means1)
            .zip(r.iter())
            .map(|((m2, m1), ri)| m2 - m1 - ri)
            .collect();

        let r_norm = dot_product(&r, &r).sqrt();
        let r2: Vec<f64> = means2.iter().zip(&means1).map(|(a, b)| a - b).collect();
        let r2_norm = dot_product(&r2, &r2).sqrt();
        let v_norm = dot_product(&v, &v).sqrt();
        let rr = dot_product(&r, &v);
        let rv_norm = rr.abs().sqrt();

        if v_norm == 0.0 {
            means0 = means1;
            break;
        }

        let alpha_s = r_norm / rv_norm;
        let alpha_s: f64 = min_step.max(alpha_s.min(max_step));

        if r_norm < OPT_TOL || r2_norm < OPT_TOL {
            if r2_norm < OPT_TOL {
                means0 = means2;
            }
            break;
        }

        // Extrapolation step
        let mut means_prime: Vec<f64> = means0
            .iter()
            .zip(&r)
            .zip(&v)
            .map(|((&m0, ri), vi)| (m0 + 2.0 * alpha_s * ri + alpha_s * alpha_s * vi).max(0.0))
            .collect();

        // Stabilization step
        if (alpha_s - 1.0).abs() > 0.01 {
            means_prime = em_update(
                &means_prime,
                te_index,
                uniq_counts,
                multi_reads,
                estimated_read_length,
            );

            if (alpha_s - max_step).abs() < 1e-10 {
                max_step = max_step0.max(max_step / m_step);
            }
        }

        if (alpha_s - max_step).abs() < 1e-10 {
            max_step *= m_step;
        }
        if min_step < 0.0 && (alpha_s - min_step).abs() < 1e-10 {
            min_step *= m_step;
        }

        means0 = means_prime;
    }

    if cur_iter >= num_iterations {
        eprintln!("EM did not converge after {} iterations", num_iterations);
    } else {
        eprintln!("EM converged at iteration {}", cur_iter);
    }

    compute_abundances(&means0, multi_reads)
}
