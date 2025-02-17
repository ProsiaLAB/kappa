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
use num_complex::Complex;
use statrs::function::gamma::ln_gamma;

use crate::utils::{spline, splint};

pub struct MieConfig {
    pub rad: f64,
    pub nangle: usize,
    pub lam: f64,
    pub nparts: usize,
    pub develop: bool,
    pub nsubr: Vec<usize>,
    pub ngaur: Vec<usize>,
    pub idis: Vec<usize>,
    pub ndis: usize,
    pub numr: usize,
    pub delta: f64,
    pub cutoff: f64,
    pub rdis: Array2<f64>,
    pub nwrdis: Array2<f64>,
    pub thmin: f64,
    pub thmax: f64,
    pub step: f64,
    pub wavelengths: Vec<f64>,
    pub cmm: Vec<Complex<f64>>,
}

pub struct MieResult {
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
    pub f_0: Array2<f64>,
}

pub enum MieError {
    InvalidSizeDistributionIndex,
    TermsOverflow,
    IntegrationAnglesOverflow,
    ScatteringAnglesOverflow,
    InvalidSplineArgument,
    DHSInvalidArgument,
    DHSInvalidDimension,
}

pub fn de_rooij_1984(miec: &mut MieConfig) -> Result<Vec<MieResult>, MieError> {
    miec.delta = 1e-8;
    miec.cutoff = 1e-8;

    miec.thmin = 180.0 * (1.0 - 0.5) / miec.nangle as f64;
    miec.thmax = 180.0 * (miec.nangle as f64 - 0.5) / miec.nangle as f64;
    miec.step = (miec.thmax - miec.thmin) / (miec.nangle as f64 - 1.0);

    miec.cmm = vec![Complex::new(0.0, 0.0); miec.nparts];
    miec.cmm[0] = Complex::new(1.0, 0.0);

    let mut rdis: Array2<f64> = Array2::zeros((miec.nparts, miec.ndis));
    rdis[[0, 0]] = miec.rad;

    let mut nwrdis: Array2<f64> = Array2::zeros((miec.nparts, miec.ndis));
    nwrdis[[0, 0]] = 1.0;

    let p1: Vec<f64> = vec![miec.rad; miec.numr];
    let p2: Vec<f64> = vec![0.0; miec.numr];
    let p3: Vec<f64> = vec![0.0; miec.numr];

    let mut f_11: Vec<f64> = vec![0.0; miec.nangle];
    let mut f_12: Vec<f64> = vec![0.0; miec.nangle];
    let mut f_33: Vec<f64> = vec![0.0; miec.nangle];
    let mut f_34: Vec<f64> = vec![0.0; miec.nangle];

    let miers = mie(&miec, &p1, &p2, &p3)?;
    for ipart in 0..miec.nparts {
        for i in 0..miec.nangle {
            f_11[i] = miers[ipart].f_0[[0, miec.nangle - 1]];
            f_12[i] = miers[ipart].f_0[[1, miec.nangle - 1]];
            f_33[i] = miers[ipart].f_0[[2, miec.nangle - 1]];
            f_34[i] = miers[ipart].f_0[[3, miec.nangle - 1]];
        }
    }
    Ok(miers)
}

fn mie(
    miec: &MieConfig,
    p1: &Vec<f64>,
    p2: &Vec<f64>,
    p3: &Vec<f64>,
) -> Result<Vec<MieResult>, MieError> {
    let mut miers: Vec<MieResult> = vec![];

    for ipart in 0..miec.nparts {
        let (rmin, rmax) = get_integration_bounds(ipart, &miec, p1, p2, p3)?;

        let m = miec.cmm[ipart].conj();

        let mier = get_scattering_matrix(&miec, ipart, m, rmin, rmax)?;

        miers.push(mier);
    }

    Ok(miers)
}

/// Find the integration bounds rmin and rmax for the integration over
/// a size distribution. These bounds are chosen such that the size
/// distribution falls below the user specified cutoff. It is essential
/// that the size distribution is normalized such that the integral
/// over all `r` is equal to one.
fn get_integration_bounds(
    ipart: usize,
    miec: &MieConfig,
    p1: &Vec<f64>,
    p2: &Vec<f64>,
    p3: &Vec<f64>,
) -> Result<(f64, f64), MieError> {
    let eps = 1e-10;

    let mut r: Vec<f64> = vec![0.0; miec.numr];

    let sef: f64;
    let ref_0: f64;
    let rref: f64;

    let rmin: f64;
    let rmax: f64;

    let idis = miec.idis[ipart];
    let p1 = p1[ipart];
    let p2 = p2[ipart];
    let p3 = p3[ipart];

    match idis {
        0 => return Ok((p1, p1)),
        1 => {
            sef = 1.0 / (p2 + 3.0).sqrt();
            ref_0 = 1.0 / (sef * sef * p2);
            rref = ref_0;
        }
        2 => {
            ref_0 = p1;
            sef = p2.sqrt();
            rref = ref_0;
        }
        3 => {
            sef = p3.sqrt();
            ref_0 = p1.max(p2) + sef;
            rref = p1.min(p2);
        }
        4 => {
            sef = (p2.ln().powi(2).exp() - 1.0).sqrt();
            ref_0 = p1 * (1.0 + sef * sef).powf(0.4);
            rref = ref_0;
        }
        5 => {
            ref_0 = p1;
            sef = ref_0.sqrt();
            rref = ref_0;
        }
        6 => {
            return Ok((p2, p3));
        }
        7 => {
            ref_0 = p2;
            sef = 2.0 * ref_0;
            rref = 0.5 * ref_0;
        }
        8 => {
            ref_0 = (p1 / (p2 * p3)).powf(p3);
            sef = 2.0 * ref_0;
            rref = 0.5 * ref_0;
        }
        9 => {
            return Ok((p1, p2));
        }
        _ => return Err(MieError::InvalidSizeDistributionIndex),
    }

    r[0] = ref_0 + sef;
    let mut r0 = ref_0;

    let mut nwithr: Vec<f64> = vec![0.0; miec.numr];

    while nwithr[0] > miec.cutoff {
        r0 = r[0];
        r[0] = 2.0 * r[0];
        get_size_distribution(
            idis,
            p1,
            p2,
            p3,
            &r,
            1,
            miec.ndis,
            &miec.rdis,
            &miec.nwrdis,
            &mut nwithr,
        )?;
    }

    let mut r1 = r[0];

    while r1 - r0 > eps {
        r[0] = 0.5 * (r0 + r1);
        get_size_distribution(
            idis,
            p1,
            p2,
            p3,
            &r,
            1,
            miec.ndis,
            &miec.rdis,
            &miec.nwrdis,
            &mut nwithr,
        )?;
        if nwithr[0] > miec.cutoff {
            r0 = r[0];
        } else {
            r1 = r[0];
        }
    }
    rmax = 0.5 * (r0 + r1);

    r1 = rref;
    r[0] = 0.0;

    while r1 > eps {
        r[0] = 0.5 * r1;
        get_size_distribution(
            idis,
            p1,
            p2,
            p3,
            &r,
            1,
            miec.ndis,
            &miec.rdis,
            &miec.nwrdis,
            &mut nwithr,
        )?;
        if nwithr[0] > miec.cutoff {
            r1 = r[0];
        } else {
            r0 = r[0];
        }
    }

    while r1 - r0 > eps {
        r[0] = 0.5 * (r0 + r1);
        get_size_distribution(
            idis,
            p1,
            p2,
            p3,
            &r,
            1,
            miec.ndis,
            &miec.rdis,
            &miec.nwrdis,
            &mut nwithr,
        )?;
        if nwithr[0] > miec.cutoff {
            r1 = r[0];
        } else {
            r0 = r[0];
        }
    }
    if r1 <= eps {
        rmin = 0.0;
    } else {
        rmin = 0.5 * (r0 + r1);
    }

    Ok((rmin, rmax))
}

fn get_size_distribution(
    idis: usize,
    p1: f64,
    p2: f64,
    p3: f64,
    r: &Vec<f64>,
    numr: usize,
    ndis: usize,
    rdis: &Array2<f64>,
    nwrdis: &Array2<f64>,
    nwithr: &mut Vec<f64>,
) -> Result<(), MieError> {
    let alpha_0: f64;
    let alpha_1: f64;

    let b0: f64;
    let b1: f64;
    let b2: f64;

    let c0: f64;

    let log_c0: f64;
    let log_c1: f64;
    let log_c2: f64;

    let flogrg: f64;
    let flogsi: f64;

    let fac: f64;
    let rg: f64;
    let aperg: f64;

    match idis {
        0 => Ok(()),
        1 => {
            // Two parameter gamms with alpha and b given
            alpha_0 = p1;
            b0 = p2;
            alpha_1 = alpha_0 + 1.0;
            log_c0 = alpha_1 * b0.ln() - ln_gamma(alpha_1);
            for i in 0..numr {
                nwithr[i] = (log_c0 + alpha_0 * r[i].ln() - b0 * r[i]).exp();
            }
            Ok(())
        }
        2 => {
            // Two parameter gamma with `p1 = reff` and `p2 = veff` given
            alpha_0 = 1.0 / p2 - 3.0;
            b0 = 1.0 / (p1 * p2);
            alpha_1 = alpha_0 + 1.0;
            log_c0 = alpha_1 * b0.ln() - ln_gamma(alpha_1);
            for i in 0..numr {
                nwithr[i] = (log_c0 + alpha_0 * r[i].ln() - b0 * r[i]).exp();
            }
            Ok(())
        }
        3 => {
            alpha_0 = 1.0 / p3 - 3.0;
            b1 = 1.0 / (p1 * p3);
            b2 = 1.0 / (p2 * p3);
            let gamlna = ln_gamma(alpha_0 + 1.0);
            log_c1 = (alpha_0 + 1.0) * b1.ln() - gamlna;
            log_c2 = (alpha_0 + 1.0) * b2.ln() - gamlna;
            for i in 0..numr {
                nwithr[i] = 0.5
                    * ((log_c1 + alpha_0 * r[i].ln() - b1 * r[i]).exp()
                        + (log_c2 + alpha_0 * r[i].ln() - b2 * r[i]).exp());
            }
            Ok(())
        }
        4 => {
            flogrg = p1.ln();
            flogsi = p2.ln().abs();
            c0 = 1.0 / ((2.0 * PI).sqrt() * flogsi);
            fac = -0.5 / (flogsi * flogsi);
            for i in 0..numr {
                nwithr[i] = c0 * (fac * (r[i].ln() - flogrg).powi(2)).exp() / r[i];
            }
            Ok(())
        }
        5 => {
            rg = p1 / (1.0 + p2).powf(2.5);
            flogrg = rg.ln();
            flogsi = (1.0 + p2).ln().sqrt();
            c0 = 1.0 / ((2.0 * PI).sqrt() * flogsi);
            fac = -0.5 / (flogsi * flogsi);
            for i in 0..numr {
                nwithr[i] = c0 * (fac * (r[i].ln() - flogrg).powi(2)).exp() / r[i];
            }
            Ok(())
        }
        6 => {
            alpha_0 = p1;
            let rmin = p2;
            let rmax = p3;
            if (alpha_0 + 1.0).abs() < 1e-10 {
                c0 = 1.0 / (rmax / rmin).ln();
            } else {
                alpha_1 = alpha_0 - 1.0;
                c0 = alpha_1 * rmax.powf(alpha_1) / ((rmax / rmin).powf(alpha_1) - 1.0);
            }
            for i in 0..numr {
                if r[i] < rmax && r[i] > rmin {
                    nwithr[i] = c0 * r[i].powf(alpha_0);
                } else {
                    nwithr[i] = 0.0;
                }
            }
            Ok(())
        }
        7 => {
            alpha_0 = p1;
            let rc = p2;
            let gamma = p3;
            b0 = alpha_0 / (gamma * rc.powf(gamma));
            aperg = (alpha_0 + 1.0) / gamma;
            log_c0 = gamma.ln() + aperg * b0.ln() - ln_gamma(aperg);
            for i in 0..numr {
                nwithr[i] = (log_c0 + alpha_0 * r[i].ln() - b0 * r[i].powf(gamma)).exp();
            }
            Ok(())
        }
        8 => {
            alpha_0 = p1;
            b0 = p2;
            let gamma = p3;
            aperg = (alpha_0 + 1.0) / gamma;
            log_c0 = gamma.ln() + aperg * b0.ln() - ln_gamma(aperg);
            for i in 0..numr {
                nwithr[i] = (log_c0 + alpha_0 * r[i].ln() - b0 * r[i].powf(gamma)).exp();
            }
            Ok(())
        }
        9 => {
            let mut rdisp: Vec<f64> = vec![0.0; ndis];
            let mut nwrdisp: Vec<f64> = vec![0.0; ndis];
            for k in 0..ndis {
                rdisp[k] = rdis[[0, k]];
                nwrdisp[k] = nwrdis[[0, k]];
            }
            let y2 = spline(&rdisp, &nwrdisp, ndis, 1e100, 1e100);
            for j in 0..numr {
                match splint(&rdisp, &nwrdisp, &y2, ndis, r[j]) {
                    Ok(n) => nwithr[j] = n,
                    Err(_) => return Err(MieError::InvalidSplineArgument),
                }
            }
            Ok(())
        }
        _ => {
            return Err(MieError::InvalidSizeDistributionIndex);
        }
    }
}

fn get_scattering_matrix(
    miec: &MieConfig,
    ipart: usize,
    m: Complex<f64>,
    rmin: f64,
    rmax: f64,
) -> Result<MieResult, MieError> {
    let mut mier = MieResult {
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

    let nfou = 0.0;
    let fac_90 = 0.0;
    let ci = Complex::new(0.0, 1.0);

    let symmetric = test_symmetry(miec.thmin, miec.thmax, miec.step);

    const RAD_FAC: f64 = PI / 180.0;

    let lambda = miec.lam;

    let rtox = 2.0 * PI / lambda;
    let fakt = lambda * lambda / (2.0 * PI);

    let nfac = 0;

    let mut w: Vec<f64> = vec![0.0; miec.numr];
    let mut r: Vec<f64> = vec![0.0; miec.numr];
    let mut nwithr: Vec<f64> = vec![0.0; miec.numr];

    let nsub: usize;
    let ngauss: usize;
    let dr: f64;

    if miec.idis[ipart] == 0 {
        w[0] = 1.0;
        r[0] = rmin;
        nwithr[0] = 1.0;

        nsub = 1;
        ngauss = 1;
        dr = 0.0;
    } else {
        dr = (rmax - rmin) / (miec.nsubr[ipart] - 1) as f64;
        // Call Gauss-Legendre quadrature routine
        // After that get the size distribution
        todo!()
    }

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
