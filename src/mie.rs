//! Mie scattering calculation
//!
//! This module contains functions to calculate Mie scattering
//! coefficients for a given set of parameters.
//!
//! # References
//! - De Rooij, W. A., & Van der Stap, C. C. (1984). The determination of the Mie

use std::f64::consts::PI;

use anyhow::{bail, Result};
use ndarray::Array2;
use num_complex::Complex;
use statrs::function::gamma::ln_gamma;

use crate::utils::{spline, splint};

pub fn de_rooij_1984(rad: f64, nangle: usize, lam: f64, f_11: Vec<f64>) {
    let develop = 0;
    let delta = 1e-8;
    let cutoff = 1e-8;

    let thmin = 180.0 * (1.0 - 0.5) / nangle as f64;
    let thmax = 180.0 * (nangle as f64 - 0.5) / nangle as f64;
    let step = (thmax - thmin) / (nangle as f64 - 1.0);
    let rdis = rad;

    mie();
}

fn get_integration_bounds(idis: usize, p1: f64, p2: f64, p3: f64) -> Result<(f64, f64)> {
    let mut r: Vec<f64> = vec![0.0; 1];

    let sef: f64;
    let ref_0: f64;
    let rref: f64;

    match idis {
        0 => return Ok((p1, p2)),
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
        _ => bail!("Illegal size distribution index"),
    }

    r[0] = ref_0 + sef;
    let r0 = ref_0;

    get_size_distribution(idis, p1, p2);

    Ok((r0, r[0]))
}

fn get_size_distribution(
    idis: usize,
    p1: f64,
    p2: f64,
    p3: f64,
    r: Vec<f64>,
    iparts: usize,
    ndis: usize,
    rdis: Vec<f64>,
    mut nwithr: Vec<f64>,
) {
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

    let ndis_max = 300;
    let ndpart = 4;
    let nwrdis: Array2<f64> = Array2::zeros((ndpart, ndis_max));
    let nwrdisp: Vec<f64> = vec![0.0; ndis_max];

    let rdis: Array2<f64> = Array2::zeros((ndpart, ndis_max));
    let rdisp: Vec<f64> = vec![0.0; ndis_max];

    let y2: Vec<f64> = vec![0.0; ndis_max];

    match idis {
        0 => {
            return;
        }
        1 => {
            // Two parameter gamms with alpha and b given
            alpha_0 = p1;
            b0 = p2;
            alpha_1 = alpha_0 + 1.0;
            log_c0 = alpha_1 * b0.ln() - ln_gamma(alpha_1);
            for i in 0..ndis {
                nwithr[i] = (log_c0 + alpha_0 * r[i].ln() - b0 * r[i]).exp();
            }
        }
        2 => {
            // Two parameter gamma with `p1 = reff` and `p2 = veff` given
            alpha_0 = 1.0 / p2 - 3.0;
            b0 = 1.0 / (p1 * p2);
            alpha_1 = alpha_0 + 1.0;
            log_c0 = alpha_1 * b0.ln() - ln_gamma(alpha_1);
            for i in 0..ndis {
                nwithr[i] = (log_c0 + alpha_0 * r[i].ln() - b0 * r[i]).exp();
            }
        }
        3 => {
            alpha_0 = 1.0 / p3 - 3.0;
            b1 = 1.0 / (p1 * p3);
            b2 = 1.0 / (p2 * p3);
            let gamlna = ln_gamma(alpha_0 + 1.0);
            log_c1 = (alpha_0 + 1.0) * b1.ln() - gamlna;
            log_c2 = (alpha_0 + 1.0) * b2.ln() - gamlna;
            for i in 0..ndis {
                nwithr[i] = 0.5
                    * ((log_c1 + alpha_0 * r[i].ln() - b1 * r[i]).exp()
                        + (log_c2 + alpha_0 * r[i].ln() - b2 * r[i]).exp());
            }
        }
        4 => {
            flogrg = p1.ln();
            flogsi = p2.ln().abs();
            c0 = 1.0 / ((2.0 * PI).sqrt() * flogsi);
            fac = -0.5 / (flogsi * flogsi);
            for i in 0..ndis {
                nwithr[i] = c0 * (fac * (r[i].ln() - flogrg).powi(2)).exp() / r[i];
            }
        }
        5 => {
            rg = p1 / (1.0 + p2).powf(2.5);
            flogrg = rg.ln();
            flogsi = (1.0 + p2).ln().sqrt();
            c0 = 1.0 / ((2.0 * PI).sqrt() * flogsi);
            fac = -0.5 / (flogsi * flogsi);
            for i in 0..ndis {
                nwithr[i] = c0 * (fac * (r[i].ln() - flogrg).powi(2)).exp() / r[i];
            }
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
            for i in 0..ndis {
                if r[i] < rmax && r[i] > rmin {
                    nwithr[i] = c0 * r[i].powf(alpha_0);
                } else {
                    nwithr[i] = 0.0;
                }
            }
        }
        7 => {
            alpha_0 = p1;
            let rc = p2;
            let gamma = p3;
            b0 = alpha_0 / (gamma * rc.powf(gamma));
            aperg = (alpha_0 + 1.0) / gamma;
            log_c0 = gamma.ln() + aperg * b0.ln() - ln_gamma(aperg);
            for i in 0..ndis {
                nwithr[i] = (log_c0 + alpha_0 * r[i].ln() - b0 * r[i].powf(gamma)).exp();
            }
        }
        8 => {
            alpha_0 = p1;
            b0 = p2;
            let gamma = p3;
            aperg = (alpha_0 + 1.0) / gamma;
            log_c0 = gamma.ln() + aperg * b0.ln() - ln_gamma(aperg);
            for i in 0..ndis {
                nwithr[i] = (log_c0 + alpha_0 * r[i].ln() - b0 * r[i].powf(gamma)).exp();
            }
        }
        9 => {
            let mut rdisp: Vec<f64> = vec![0.0; ndis];
            let mut nwrdisp: Vec<f64> = vec![0.0; ndis];
            for k in 0..ndis {
                rdisp[k] = rdis[[0, k]];
                nwrdisp[k] = nwrdis[[0, k]];
            }
            spline(xv, yv, n, yp1, ypn);
            for j in 0..numr {
                splint(xv, yv, y2, x, jl);
            }
        }
        _ => {
            return;
        }
    }
}

fn mie(wavel: f64, real: f64, imag: f64) {
    let lambda = wavel;
    let nr = real;
    let ni = imag;
    let m = Complex::new(nr, -ni);

    get_scattering_matrix();
}

fn get_scattering_matrix(
    m: Complex<f64>,
    lambda: f64,
    p1: f64,
    p2: f64,
    p3: f64,
    delta: f64,
    thmin: f64,
    thmax: f64,
    step: f64,
    mut nangle: usize,
) -> Result<(Vec<f64>, Vec<f64>, Array2<f64>, Vec<f64>, usize)> {
    let mut f: Array2<f64> = Array2::zeros((3, 3));

    let c_sca = 0.0;
    let c_ext = 0.0;
    let numpar = 0.0;
    let g: f64 = 0.0;
    let reff: f64 = 0.0;
    let nfou: f64 = 0.0;
    let fac_90: f64 = 0.0;
    let ci: Complex<f64> = Complex::new(0.0, 1.0);

    let symmetric = test_symmetry(thmin, thmax, step);

    const RAD_FAC: f64 = PI / 180.0;

    let rtox = 2.0 * PI / lambda;
    let fakt = lambda * lambda / (2.0 * PI);

    let nfac = 0;

    let w = 1.0;
    let r = p1;
    let nwithr = 1.0;

    let mut u: Vec<f64> = vec![0.0; nangle];
    let mut wth: Vec<f64> = vec![0.0; nangle];
    let mut miec: Vec<f64> = vec![0.0; nangle];

    Ok((u, wth, f, miec, nangle))
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
