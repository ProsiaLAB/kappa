//! Mie scattering calculation
//!
//! This module contains functions to calculate Mie scattering
//! coefficients for a given set of parameters.
//!
//! # References
//! - De Rooij, W. A., & Van der Stap, C. C. (1984).

use std::f64::consts::PI;

use anyhow::anyhow;
use anyhow::Result;
use num_complex::{Complex, Complex64, ComplexFloat};

use crate::types::{CVector, RMatrix, RVector};

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

pub fn de_rooij_1984(miec: &MieConfig) -> Result<MieResult> {
    let m = miec.cmm.conj();
    let mier = get_scattering_matrix(miec, m)?;
    Ok(mier)
}

fn get_scattering_matrix(miec: &MieConfig, m: Complex64) -> Result<MieResult> {
    let mut mier = MieResult {
        u: RVector::zeros(miec.nangle),
        wth: RVector::zeros(miec.nangle),
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
        f_0: RMatrix::zeros((4, miec.nangle)),
        f_11: RVector::zeros(miec.nangle),
        f_12: RVector::zeros(miec.nangle),
        f_33: RVector::zeros(miec.nangle),
        f_34: RVector::zeros(miec.nangle),
    };

    let mut fac_90 = 0.0;
    let ci = Complex::new(0.0, 1.0);

    let symmetric = test_symmetry(miec.thmin, miec.thmax, miec.step);

    const RAD_FAC: f64 = PI / 180.0;

    let lambda = miec.lam;

    let rtox = 2.0 * PI / lambda;
    let fakt = lambda * lambda / (2.0 * PI);

    let nfac = 0;

    let w = 1.0;
    let r = miec.rad;
    let nwithr = 1.0;

    let sw = nwithr * w;
    let x = rtox * r;
    let nmax = (x + 4.05 * x.powf(1.0 / 3.0) + 2.0) as usize;
    let nfi = nmax + 60;
    let zabs = x * m.abs();
    let nd = (zabs + 4.05 * zabs.powf(1.0 / 3.0) + 70.0) as usize;

    let size = nd.max(nfi).max(nmax);

    let mut facf = RVector::zeros(size);
    let mut facb = RVector::zeros(size);

    let (fi, chi, d) = fichid(m, x, nfi, nmax, nd);
    let (an, bn) = anbn(m, x, fi, chi, d, nmax);

    if nmax > nfac {
        for n in nfac + 1..nmax {
            let nf = n as f64;
            facf[n] = (2.0 * nf + 1.0) / (nf * (nf + 1.0));
            facb[n] = facf[n];
            if n % 2 == 1 {
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
        aux = (2.0 * nf + 1.0) * (an[n].norm() + bn[n].norm()).abs();
        c_sca_sum += aux;
        c_ext_sum += (2.0 * nf + 1.0) * (an[n] + bn[n]).re;
        if aux < miec.delta {
            nstop = n;
            break;
        }
    }

    let nfou = nstop;

    if nfou > nmax {
        eprintln!("WARNING: `mie` sum not converged for scattering cross-section");
        eprintln!("         with radius = {} and size parameter = {}", r, x);
        eprintln!("         size distribution `nr` = {}", 0);
        eprintln!("         Re(m) = {} and Im(m) = {}", m.re, m.im);
        eprintln!(
            "         apriori estimate of number of `mie` terms = {}",
            nmax
        );
        eprintln!("         Term {} for `c_sca` was {}", nmax, aux);
        eprintln!(
            "         should have been less than `delta` = {}",
            miec.delta
        );
        eprintln!("         the apriori estimate will be used instead.");
    }

    let nangle = ((miec.thmax - miec.thmin) / miec.step) as usize + 1;
    if nangle > 6000 {
        return Err(anyhow!("ScatteringAnglesOverflow"));
    }
    mier.u = RVector::zeros(nangle);
    mier.wth = RVector::zeros(nangle);
    let wfac = 2.0 / nangle as f64;
    for iang in 0..nangle {
        let iangf = iang as f64;
        let th = miec.thmin + iangf * miec.step;
        mier.u[nangle - 1 - iang] = (RAD_FAC * th).cos();
        mier.wth[iang] = wfac;
    }

    mier.numpar += sw;
    mier.g += sw * r.powi(3);

    let mut nhalf: usize = 0;
    if symmetric {
        if nangle % 2 == 1 {
            nhalf = nangle.div_ceil(2);
            fac_90 = 0.5;
        } else {
            nhalf = nangle / 2;
        }

        for j in 0..nhalf {
            let (pi, tau) = pitau(mier.u[j], nmax);
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
            mier.f_0[[0, j]] += sw * (s_plus_f * conj_s_plus_f + s_min_f * conj_s_min_f).re;
            mier.f_0[[1, j]] -= sw * (s_min_f * conj_s_plus_f + s_plus_f * conj_s_min_f).re;
            mier.f_0[[2, j]] += sw * (s_plus_f * conj_s_plus_f - s_min_f * conj_s_min_f).re;
            mier.f_0[[3, j]] += (ci * sw * (s_min_f * conj_s_plus_f - s_plus_f * conj_s_min_f)).re;

            // Backward scattering elements
            mier.f_0[[0, k]] += sw * (s_plus_b * conj_s_plus_b + s_min_b * conj_s_min_b).re;
            mier.f_0[[1, k]] -= sw * (s_min_b * conj_s_plus_b + s_plus_b * conj_s_min_b).re;
            mier.f_0[[2, k]] += sw * (s_plus_b * conj_s_plus_b - s_min_b * conj_s_min_b).re;
            mier.f_0[[3, k]] += (ci * sw * (s_min_b * conj_s_plus_b - s_plus_b * conj_s_min_b)).re;
        }
    } else {
        for j in 0..nangle {
            let (pi, tau) = pitau(mier.u[j], nmax);
            let mut s_plus_f = Complex::new(0.0, 0.0);
            let mut s_min_f = Complex::new(0.0, 0.0);

            for n in 0..nfou {
                s_plus_f += facf[n] * (an[n] + bn[n]) * (pi[n] + tau[n]);
                s_min_f += facf[n] * (an[n] - bn[n]) * (pi[n] - tau[n]);
            }
            let conj_s_plus_f = s_plus_f.conj();
            let conj_s_min_f = s_min_f.conj();

            // Forward scattering elements
            mier.f_0[[0, j]] += sw * (s_plus_f * conj_s_plus_f + s_min_f * conj_s_min_f).re;
            mier.f_0[[1, j]] -= sw * (s_min_f * conj_s_plus_f + s_plus_f * conj_s_min_f).re;
            mier.f_0[[2, j]] += sw * (s_plus_f * conj_s_plus_f - s_min_f * conj_s_min_f).re;
            mier.f_0[[3, j]] += (ci * sw * (s_min_f * conj_s_plus_f - s_plus_f * conj_s_min_f)).re;
        }
    }

    mier.c_sca += sw * c_sca_sum;
    mier.c_ext += sw * c_ext_sum;

    for j in 0..nangle {
        for k in 0..4 {
            mier.f_0[[k, j]] /= 2.0 * mier.c_sca;
        }
    }

    if symmetric {
        for k in 0..4 {
            mier.f_0[[k, nhalf - 1]] *= fac_90;
        }
    }

    mier.g *= PI;
    mier.c_sca *= fakt;
    mier.c_ext *= fakt;
    mier.q_sca = mier.c_sca / mier.g;
    mier.q_ext = mier.c_ext / mier.g;
    mier.albedo = mier.c_sca / mier.c_ext;
    mier.volume = (4.0 / 3.0) * PI * mier.reff;
    mier.reff = PI * mier.reff / mier.g;
    mier.xeff = rtox * mier.reff;

    for i in 0..miec.nangle {
        mier.f_11[i] = mier.f_0[[0, miec.nangle - 1]];
        mier.f_12[i] = mier.f_0[[1, miec.nangle - 1]];
        mier.f_33[i] = mier.f_0[[2, miec.nangle - 1]];
        mier.f_34[i] = mier.f_0[[3, miec.nangle - 1]];
    }

    Ok(mier)
}

fn test_symmetry(thmin: f64, thmax: f64, step: f64) -> bool {
    let eps = 1e-6;
    let heps = 0.5 * eps;

    ((180.0 - thmin - thmax).abs() < eps) && ((thmax - thmin + heps).rem_euclid(step) < eps)
}

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
        psi[n] = 1.0 / (2.0 * nf + 1.0) / x - psi[n + 1];
    }

    for n in (0..nd).rev() {
        let nf = n as f64;
        let zn_1 = (nf + 1.0) * perz;
        d[n] = zn_1 - 1.0 / (d[n] + zn_1);
    }

    psi[0] = sinx;
    let psi_1 = psi[0] * perx - cosx;
    if psi_1.abs() > 1e-4 {
        psi[1] = psi_1;
        for n in 2..nmax {
            psi[n] *= psi[n - 1];
        }
    } else {
        for n in 1..nmax {
            psi[n] *= psi[n - 1];
        }
    }

    chi[0] = cosx;
    chi[1] = chi[0] * perx + sinx;
    for n in 1..nmax - 1 {
        let nf = n as f64;
        chi[n + 1] = (2.0 * nf + 1.0) * chi[n] * perx - chi[n - 1];
    }

    (psi, chi, d)
}

fn anbn(
    m: Complex64,
    x: f64,
    psi: RVector,
    chi: RVector,
    d: CVector,
    nmax: usize,
) -> (CVector, CVector) {
    let perm = 1.0 / m;
    let perx = 1.0 / x;

    let mut an = CVector::zeros(nmax);
    let mut bn = CVector::zeros(nmax);

    for n in 0..nmax {
        let nf = n as f64;
        let zn = Complex::new(psi[n], chi[n]);
        let znm_1 = Complex::new(psi[n], chi[n]);
        let xn = nf * perx;
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

    for n in 2..(nmax - 1) {
        let nf = n as f64;
        pi[n + 1] = (nf + 1.0) / (nf) * delta + u * pi[n];
        delta = u * pi[n + 1] - pi[n];
        tau[n + 1] = (nf + 1.0) * delta - pi[n];
    }

    (pi, tau)
}
