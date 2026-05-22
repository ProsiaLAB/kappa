//! Mie scattering calculation
//!
//! This module contains functions to calculate Mie scattering
//! coefficients for a given set of parameters.
//!
//! # References
//! - De Rooij, W. A., & Van der Stap, C. C. (1984).

use std::f64::consts::PI;

use anyhow::Result;
use anyhow::anyhow;
use extensions::types::{CVector, RMatrix, RVector};
use num_complex::{Complex, Complex64, ComplexFloat};

#[derive(Debug)]
pub struct MieConfig {
    pub nangle: usize,
    pub delta: f64,
    pub thmin: f64,
    pub thmax: f64,
    pub step: f64,
    pub lam: f64,
    pub cmm: Complex64,
    pub rad: f64,
}

#[derive(Debug)]
pub struct MieResult {
    pub u: RVector,
    pub wth: RVector,
    pub c_sca: f64,
    pub c_ext: f64,
    pub q_sca: f64,
    pub q_ext: f64,
    pub albedo: f64,
    pub g: f64,
    pub reff: f64,
    pub xeff: f64,
    pub numpar: f64,
    pub volume: f64,
    pub f_0: RMatrix,
    pub f_11: RVector,
    pub f_12: RVector,
    pub f_33: RVector,
    pub f_34: RVector,
}

// pub enum MieError {
//     ScatteringAnglesOverflow,
// }

/// Calculate Mie scattering coefficients using the method of De Rooij & Van der Stap (1984).
///
/// # Errors
/// - `ScatteringAnglesOverflow`: The number of scattering angles exceeds the maximum allowed (6000).
/// - `MieSumNotConverged`: The Mie sum did not converge for the given parameters, and the apriori estimate of the number of terms will be used instead.
/// - `InvalidMieConfig`: The provided Mie configuration is invalid (e.g., negative radius, non-positive wavelength).
pub fn de_rooij_1984(mie_cfg: &MieConfig) -> Result<MieResult> {
    let m = mie_cfg.cmm.conj();
    let mie_result = get_scattering_matrix(mie_cfg, m)?;
    Ok(mie_result)
}

#[allow(clippy::similar_names)]
fn get_scattering_matrix(mie_cfg: &MieConfig, m: Complex64) -> Result<MieResult> {
    const RAD_FAC: f64 = PI / 180.0;

    let mut mie_result = MieResult {
        u: RVector::zeros(mie_cfg.nangle),
        wth: RVector::zeros(mie_cfg.nangle),
        c_sca: 0.0,
        c_ext: 0.0,
        q_sca: 0.0,
        q_ext: 0.0,
        albedo: 0.0,
        g: 0.0,
        reff: 0.0,
        xeff: 0.0,
        numpar: 0.0,
        volume: 0.0,
        f_0: RMatrix::zeros((4, mie_cfg.nangle)),
        f_11: RVector::zeros(mie_cfg.nangle),
        f_12: RVector::zeros(mie_cfg.nangle),
        f_33: RVector::zeros(mie_cfg.nangle),
        f_34: RVector::zeros(mie_cfg.nangle),
    };

    let mut fac_90 = 1.0;
    let ci = Complex::new(0.0, 1.0);

    let symmetric = test_symmetry(mie_cfg.thmin, mie_cfg.thmax, mie_cfg.step);

    let lambda = mie_cfg.lam;

    let rtox = 2.0 * PI / lambda;
    let fakt = lambda * lambda / (2.0 * PI);

    let nfac = 0;

    let w = 1.0;
    let r = mie_cfg.rad;
    let nwithr = 1.0;

    let sw = nwithr * w;
    let x = rtox * r;
    let nmax = (x + 4.05 * x.powf(1.0 / 3.0) + 20.0) as usize;
    let nfi = nmax + 60;
    let zabs = x * m.abs();
    let nd = (zabs + 4.05 * zabs.powf(1.0 / 3.0) + 70.0) as usize;

    let size = nd.max(nfi).max(nmax);

    let mut facf = RVector::zeros(size);
    let mut facb = RVector::zeros(size);

    let (fi, chi, d) = fichid(m, x, nfi, nmax, nd);
    let (an, bn) = anbn(m, x, &fi, &chi, &d, nmax);

    if nmax > nfac {
        for n in nfac..nmax {
            let nf = n as f64 + 1.0;
            facf[n] = (2.0 * nf + 1.0) / (nf * (nf + 1.0));
            facb[n] = facf[n];
            if n % 2 != 1 {
                facb[n] = -facb[n];
            }
        }
    }

    let mut c_ext_sum = 0.0;
    let mut c_sca_sum = 0.0;
    let mut nstop = nmax;
    let mut aux: f64 = 0.0;

    for n in 0..nmax {
        let nf = n as f64;
        let afac = an[n] * an[n].conj();
        let bfac = bn[n] * bn[n].conj();
        aux = (2.0 * (nf + 1.0) + 1.0) * (afac + bfac).abs();
        c_sca_sum += aux;
        c_ext_sum += (2.0 * (nf + 1.0) + 1.0) * (an[n] + bn[n]).re;
        if aux < mie_cfg.delta {
            nstop = n;
            break;
        }
    }

    let nfou = nstop + 1;

    if nfou > nmax {
        eprintln!("WARNING: `mie` sum not converged for scattering cross-section");
        eprintln!("         with radius = {r} and size parameter = {x}");
        eprintln!("         size distribution `nr` = {}", 0);
        eprintln!("         Re(m) = {} and Im(m) = {}", m.re, m.im);
        eprintln!("         apriori estimate of number of `mie` terms = {nmax}");
        eprintln!("         Term {nmax} for `c_sca` was {aux}");
        eprintln!(
            "         should have been less than `delta` = {}",
            mie_cfg.delta
        );
        eprintln!("         the apriori estimate will be used instead.");
    }

    let nangle = ((mie_cfg.thmax - mie_cfg.thmin) / mie_cfg.step) as usize + 1;
    if nangle > 6000 {
        return Err(anyhow!("ScatteringAnglesOverflow"));
    }
    mie_result.u = RVector::zeros(nangle);
    mie_result.wth = RVector::zeros(nangle);
    let wfac = 2.0 / nangle as f64;
    for iang in 0..nangle {
        let iangf = iang as f64;
        let th = mie_cfg.thmin + iangf * mie_cfg.step;
        mie_result.u[nangle - 1 - iang] = (RAD_FAC * th).cos();
        mie_result.wth[iang] = wfac;
    }

    mie_result.numpar += sw;
    mie_result.g += sw * r.powi(2);
    mie_result.reff += sw * r.powi(3);

    let mut nhalf: usize = 0;
    if symmetric {
        if nangle % 2 == 1 {
            nhalf = nangle.div_ceil(2);
            fac_90 = 0.5;
        } else {
            nhalf = nangle / 2;
        }

        for j in 0..nhalf {
            let (pi, tau) = pitau(mie_result.u[j], nmax);
            let mut s_plus_f = Complex::new(0.0, 0.0);
            let mut s_min_f = Complex::new(0.0, 0.0);
            let mut s_plus_b = Complex::new(0.0, 0.0);
            let mut s_min_b = Complex::new(0.0, 0.0);

            for n in 0..nfou {
                s_plus_f += facf[n] * (an[n] + bn[n]) * (pi[n] + tau[n]);
                s_min_f += facf[n] * (an[n] - bn[n]) * (pi[n] - tau[n]);
                s_plus_b += facb[n] * (an[n] + bn[n]) * (pi[n] - tau[n]);
                s_min_b += facb[n] * (an[n] - bn[n]) * (pi[n] + tau[n]);
            }
            let conj_s_plus_f = s_plus_f.conj();
            let conj_s_min_f = s_min_f.conj();
            let conj_s_plus_b = s_plus_b.conj();
            let conj_s_min_b = s_min_b.conj();

            let k = nangle - 1 - j;

            // Forward scattering elements
            mie_result.f_0[[0, j]] += sw * (s_plus_f * conj_s_plus_f + s_min_f * conj_s_min_f).re;
            mie_result.f_0[[1, j]] -= sw * (s_min_f * conj_s_plus_f + s_plus_f * conj_s_min_f).re;
            mie_result.f_0[[2, j]] += sw * (s_plus_f * conj_s_plus_f - s_min_f * conj_s_min_f).re;
            mie_result.f_0[[3, j]] +=
                (ci * sw * (s_min_f * conj_s_plus_f - s_plus_f * conj_s_min_f)).re;

            // Backward scattering elements
            mie_result.f_0[[0, k]] += sw * (s_plus_b * conj_s_plus_b + s_min_b * conj_s_min_b).re;
            mie_result.f_0[[1, k]] -= sw * (s_min_b * conj_s_plus_b + s_plus_b * conj_s_min_b).re;
            mie_result.f_0[[2, k]] += sw * (s_plus_b * conj_s_plus_b - s_min_b * conj_s_min_b).re;
            mie_result.f_0[[3, k]] +=
                (ci * sw * (s_min_b * conj_s_plus_b - s_plus_b * conj_s_min_b)).re;
        }
    } else {
        for j in 0..nangle {
            let (pi, tau) = pitau(mie_result.u[j], nmax);
            let mut s_plus_f = Complex::new(0.0, 0.0);
            let mut s_min_f = Complex::new(0.0, 0.0);

            for n in 0..nfou {
                s_plus_f += facf[n] * (an[n] + bn[n]) * (pi[n] + tau[n]);
                s_min_f += facf[n] * (an[n] - bn[n]) * (pi[n] - tau[n]);
            }
            let conj_s_plus_f = s_plus_f.conj();
            let conj_s_min_f = s_min_f.conj();

            // Forward scattering elements
            mie_result.f_0[[0, j]] += sw * (s_plus_f * conj_s_plus_f + s_min_f * conj_s_min_f).re;
            mie_result.f_0[[1, j]] -= sw * (s_min_f * conj_s_plus_f + s_plus_f * conj_s_min_f).re;
            mie_result.f_0[[2, j]] += sw * (s_plus_f * conj_s_plus_f - s_min_f * conj_s_min_f).re;
            mie_result.f_0[[3, j]] +=
                (ci * sw * (s_min_f * conj_s_plus_f - s_plus_f * conj_s_min_f)).re;
        }
    }

    mie_result.c_sca += sw * c_sca_sum;
    mie_result.c_ext += sw * c_ext_sum;

    for j in 0..nangle {
        for k in 0..4 {
            mie_result.f_0[[k, j]] /= 2.0 * mie_result.c_sca;
        }
    }

    if symmetric {
        for k in 0..4 {
            mie_result.f_0[[k, nhalf - 1]] *= fac_90;
        }
    }

    mie_result.g *= PI;
    mie_result.c_sca *= fakt;
    mie_result.c_ext *= fakt;
    mie_result.q_sca = mie_result.c_sca / mie_result.g;
    mie_result.q_ext = mie_result.c_ext / mie_result.g;
    mie_result.albedo = mie_result.c_sca / mie_result.c_ext;
    mie_result.volume = (4.0 / 3.0) * PI * mie_result.reff;
    mie_result.reff = PI * mie_result.reff / mie_result.g;
    mie_result.xeff = rtox * mie_result.reff;

    for i in 0..mie_cfg.nangle {
        mie_result.f_11[i] = mie_result.f_0[[0, mie_cfg.nangle - i - 1]];
        mie_result.f_12[i] = mie_result.f_0[[1, mie_cfg.nangle - i - 1]];
        mie_result.f_33[i] = mie_result.f_0[[2, mie_cfg.nangle - i - 1]];
        mie_result.f_34[i] = mie_result.f_0[[3, mie_cfg.nangle - i - 1]];
    }

    Ok(mie_result)
}

fn test_symmetry(thmin: f64, thmax: f64, step: f64) -> bool {
    let eps = 1e-6;
    let heps = 0.5 * eps;

    ((180.0 - thmin - thmax).abs() < eps) && ((thmax - thmin + heps).rem_euclid(step) < eps)
}

#[allow(clippy::similar_names)]
fn fichid(
    m: Complex64,
    x: f64,
    nchi: usize,
    nmax: usize,
    nd: usize,
) -> (RVector, RVector, CVector) {
    let z = m * x;
    let perz = 1.0 / z;
    let perx = 1.0 / x;

    let sinx = x.sin();
    let cosx = x.cos();

    let mut psi = RVector::zeros(nchi + 1);
    let mut chi = RVector::zeros(nmax + 2);
    let mut d = CVector::zeros(nd);

    for n in (0..nchi).rev() {
        let nf = n as f64;
        psi[n] = 1.0 / ((2.0 * nf + 1.0) / x - psi[n + 1]);
    }

    for n in (0..nd - 1).rev() {
        let nf = n as f64;
        let zn_1 = (nf + 2.0) * perz;
        d[n] = zn_1 - 1.0 / (d[n + 1] + zn_1);
    }

    psi[0] = sinx;
    let psi_1 = psi[0] * perx - cosx;
    if psi_1.abs() > 1e-4 {
        psi[1] = psi_1;
        for n in 2..=nmax {
            psi[n] *= psi[n - 1];
        }
    } else {
        for n in 1..=nmax {
            psi[n] *= psi[n - 1];
        }
    }

    chi[0] = cosx;
    chi[1] = chi[0] * perx + sinx;
    for n in 1..nmax {
        let nf = n as f64;
        chi[n + 1] = (2.0 * nf + 1.0) * chi[n] * perx - chi[n - 1];
    }

    (psi, chi, d)
}

#[allow(clippy::similar_names)]
fn anbn(
    m: Complex64,
    x: f64,
    psi: &RVector,
    chi: &RVector,
    d: &CVector,
    nmax: usize,
) -> (CVector, CVector) {
    let perm = 1.0 / m;
    let perx = 1.0 / x;

    let mut an = CVector::zeros(nmax);
    let mut bn = CVector::zeros(nmax);

    for n in 0..nmax {
        let nf = n as f64;
        let zn = Complex::new(psi[n + 1], chi[n + 1]);
        let znm_1 = Complex::new(psi[n], chi[n]);
        let xn = (nf + 1.0) * perx;
        let save_a = d[n] * perm + xn;
        an[n] = (save_a * psi[n + 1] - psi[n]) / (save_a * zn - znm_1);
        let save_b = d[n] * m + xn;
        bn[n] = (save_b * psi[n + 1] - psi[n]) / (save_b * zn - znm_1);
    }

    (an, bn)
}

fn pitau(u: f64, nmax: usize) -> (RVector, RVector) {
    let mut pi = RVector::zeros(nmax);
    let mut tau = RVector::zeros(nmax);

    pi[0] = 1.0;
    pi[1] = 3.0 * u;
    let mut delta = 3.0 * u * u - 1.0;

    tau[0] = u;
    tau[1] = 2.0 * delta - 1.0;

    for n in 1..(nmax - 1) {
        let nf = n as f64;
        pi[n + 1] = (nf + 2.0) / (nf + 1.0) * delta + u * pi[n];
        delta = u * pi[n + 1] - pi[n];
        tau[n + 1] = (nf + 2.0) * delta - pi[n];
    }

    (pi, tau)
}
