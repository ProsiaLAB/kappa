//! Mie scattering calculation
//!
//! This module contains functions to calculate Mie scattering
//! coefficients for a given set of parameters.
//!
//! # References
//! - De Rooij, W. A., & Van der Stap, C. C. (1984).

use std::f64::consts::PI;

use anyhow::Result;
use ndarray::Array2;
use num_complex::{Complex, ComplexFloat};

pub struct MieConfig {
    nangle: usize,
    delta: f64,
    thmin: f64,
    thmax: f64,
    step: f64,
    lam: f64,
    cmm: Complex<f64>,
    rad: f64,
}

pub struct MieResult {
    u: Vec<f64>,
    wth: Vec<f64>,
    c_sca: f64,
    c_ext: f64,
    q_sca: f64,
    q_ext: f64,
    albedo: f64,
    g: f64,
    reff: f64,
    xeff: f64,
    numpar: f64,
    volume: f64,
    f_0: Array2<f64>,
}

pub enum MieError {
    ScatteringAnglesOverflow,
}

pub fn de_rooij_1984(miec: &MieConfig) -> Result<MieResult, MieError> {
    // miec.delta = 1e-8;
    // miec.cutoff = 1e-8;

    // miec.thmin = 180.0 * (1.0 - 0.5) / miec.nangle as f64;
    // miec.thmax = 180.0 * (miec.nangle as f64 - 0.5) / miec.nangle as f64;
    // miec.step = (miec.thmax - miec.thmin) / (miec.nangle as f64 - 1.0);

    // miec.cmm = vec![Complex::new(0.0, 0.0); miec.nparts];
    // miec.cmm[0] = Complex::new(1.0, 0.0);

    // let mut rdis: Array2<f64> = Array2::zeros((miec.nparts, miec.ndis));
    // rdis[[0, 0]] = miec.rad;

    // let mut nwrdis: Array2<f64> = Array2::zeros((miec.nparts, miec.ndis));
    // nwrdis[[0, 0]] = 1.0;

    let mut f_11: Vec<f64> = vec![0.0; miec.nangle];
    let mut f_12: Vec<f64> = vec![0.0; miec.nangle];
    let mut f_33: Vec<f64> = vec![0.0; miec.nangle];
    let mut f_34: Vec<f64> = vec![0.0; miec.nangle];

    let mier = mie(&miec)?;

    for i in 0..miec.nangle {
        f_11[i] = mier.f_0[[0, miec.nangle - 1]];
        f_12[i] = mier.f_0[[1, miec.nangle - 1]];
        f_33[i] = mier.f_0[[2, miec.nangle - 1]];
        f_34[i] = mier.f_0[[3, miec.nangle - 1]];
    }

    Ok(mier)
}

fn mie(miec: &MieConfig) -> Result<MieResult, MieError> {
    let m = miec.cmm.conj();

    let mier = get_scattering_matrix(&miec, m)?;

    Ok(mier)
}

fn get_scattering_matrix(miec: &MieConfig, m: Complex<f64>) -> Result<MieResult, MieError> {
    let mut mier = MieResult {
        u: vec![],
        wth: vec![],
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
        f_0: Array2::zeros((3, 3)),
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
    let nfi = (nmax + 60) as usize;
    let zabs = x * m.abs();
    let nd = (zabs + 4.05 * zabs.powf(1.0 / 3.0) + 70.0) as usize;

    let size = nd.max(nfi).max(nmax);

    let mut facf: Vec<f64> = vec![0.0; size];
    let mut facb: Vec<f64> = vec![0.0; size];

    let (fi, chi, d) = fichid(m, x, nfi, nmax, nd);
    let (an, bn) = anbn(m, x, fi, chi, d, nmax);

    if nmax > nfac {
        for n in nfac + 1..nmax {
            facf[n] = (2.0 * n as f64 + 1.0) / (n as f64 * (n as f64 + 1.0));
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
        aux = (2.0 * n as f64 + 1.0) * (an[n].norm() + bn[n].norm()).abs();
        c_sca_sum += aux;
        c_ext_sum += (2.0 * n as f64 + 1.0) * (an[n] + bn[n]).re;
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
        return Err(MieError::ScatteringAnglesOverflow);
    }
    mier.u = vec![0.0; nangle];
    mier.wth = vec![0.0; nangle];
    let wfac = 2.0 / nangle as f64;
    for iang in 0..nangle {
        let th = miec.thmin + iang as f64 * miec.step;
        mier.u[nangle - 1 - iang] = (RAD_FAC * th).cos();
        mier.wth[iang] = wfac;
    }

    mier.numpar += sw;
    mier.g += sw * r.powi(3);

    let mut nhalf: usize = 0;
    if symmetric {
        if nangle % 2 == 1 {
            nhalf = (nangle + 1) / 2;
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

    Ok(mier)
}

fn test_symmetry(thmin: f64, thmax: f64, step: f64) -> bool {
    let eps = 1e-6;
    let heps = 0.5 * eps;

    if ((180.0 - thmin - thmax).abs() < eps) && ((thmax - thmin + heps).rem_euclid(step) < eps) {
        true
    } else {
        false
    }
}

fn fichid(
    m: Complex<f64>,
    x: f64,
    nchi: usize,
    nmax: usize,
    nd: usize,
) -> (Vec<f64>, Vec<f64>, Vec<Complex<f64>>) {
    let z = m * x;
    let perz = 1.0 / z;
    let perx = 1.0 / x;

    let sinx = x.sin();
    let cosx = x.cos();

    let mut psi: Vec<f64> = vec![0.0; nchi + 1];
    let mut chi: Vec<f64> = vec![0.0; nmax + 2];
    let mut d: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nd];

    for n in (0..nchi).rev() {
        psi[n] = 1.0 / (2.0 * n as f64 + 1.0) / x - psi[n + 1];
    }

    for n in (0..nd).rev() {
        let zn_1 = (n as f64 + 1.0) * perz;
        d[n] = zn_1 - 1.0 / (d[n + 1] + zn_1);
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
        chi[n + 1] = (2.0 * n as f64 + 1.0) * chi[n] * perx - chi[n - 1];
    }

    return (psi, chi, d);
}

fn anbn(
    m: Complex<f64>,
    x: f64,
    psi: Vec<f64>,
    chi: Vec<f64>,
    d: Vec<Complex<f64>>,
    nmax: usize,
) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
    let perm = 1.0 / m;
    let perx = 1.0 / x;

    let mut an: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nmax];
    let mut bn: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nmax];

    for n in 0..nmax {
        let zn = Complex::new(psi[n], chi[n]);
        let znm_1 = Complex::new(psi[n - 1], chi[n - 1]);
        let xn = n as f64 * perx;
        let save_a = d[n] * perm + xn;
        an[n] = (save_a * psi[n] - psi[n - 1]) / (save_a * zn - znm_1);
        let save_b = d[n] * m + xn;
        bn[n] = (save_b * psi[n] - psi[n - 1]) / (save_b * zn - znm_1);
    }

    return (an, bn);
}

fn pitau(u: f64, nmax: usize) -> (Vec<f64>, Vec<f64>) {
    let mut pi: Vec<f64> = vec![0.0; nmax];
    let mut tau: Vec<f64> = vec![0.0; nmax];

    pi[0] = 1.0;
    pi[1] = 3.0 * u;
    let mut delta = 3.0 * u * u - 1.0;

    tau[0] = u;
    tau[1] = 2.0 * delta - 1.0;

    for n in 1..nmax {
        pi[n + 1] = (n as f64 + 1.0) / (n as f64) * delta + u * pi[n];
        delta = u * pi[n + 1] - pi[n];
        tau[n + 1] = (n as f64 + 1.0) * delta - pi[n];
    }

    return (pi, tau);
}
