//! Computes light scattering properties of randomly oriented
//! fractal dust aggregates by means of the modified mean field theory developed
//! in Tazaki & Tanaka (2018). This code is also capable of computing the light
//! scattering solution based on the Rayleigh-Gans-Debye theory and the Mean field
//! theory.

use anyhow::bail;
use anyhow::Result;
use num_complex::Complex;

pub fn mean_scattering() {
    // maxwell_garnett_mixing();
    structure_factor_integration();
    lorenz_mie();
    todo!()
}

/// This function finds effective refractive index of an aggregates
/// by using Maxwell-Garnett mixing rule.
///
/// * eps_1 : Dielectric constant of the first material.
/// * eps_2 : Dielectric constant of vacuum.
/// * f1 : Volume filling factor of the first material.
/// * f2 : Volume filling factor of the second material (porosity).
fn maxwell_garnett_mixing(refrel: Complex<f64>, f1: f64) -> Complex<f64> {
    let eps_1 = refrel * refrel;
    let eps_2 = Complex::new(1.0, 0.0);
    let mg = eps_2 * (2.0 * f1 * (eps_1 - eps_2) + eps_1 + 2.0 * eps_2)
        / (eps_1 + 2.0 * eps_2 - f1 * (eps_1 - eps_2));

    mg.sqrt()
}

fn structure_factor_integration(d: f64, q: f64, r_g: f64) -> f64 {
    let integration_grid_size = 100_000;
    let normalization_factor = 0.5; // Botet+ 1997

    // Upper limit of integration
    let eta: f64 = 25.0;
    let xmax = (eta / normalization_factor).powf(1.0 / d);

    let mut s_q: f64;
    if q == 0.0 {
        s_q = 1.0;
    } else {
        s_q = 0.0;
        let h = xmax / (2.0 * integration_grid_size as f64);
        for i in 0..integration_grid_size - 2 {
            let x0 = 2.0 * i as f64 * h;
            let x1 = (2.0 * i as f64 + 1.0) * h;
            let x2 = (2.0 * i as f64 + 2.0) * h;
            s_q += h
                * (structure_factor_fn(x0, normalization_factor, d, q, r_g)
                    + 4.0 * structure_factor_fn(x1, normalization_factor, d, q, r_g)
                    + structure_factor_fn(x2, normalization_factor, d, q, r_g))
                / 3.0;
        }
        s_q *= normalization_factor * d / q / r_g;
    }
    if s_q < 0.0 {
        s_q = 0.0;
    }

    s_q
}

fn structure_factor_fn(x: f64, c: f64, d: f64, q: f64, r_g: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else {
        let val = (q * r_g * x).sin() * (-c * c.powf(d)).exp() * x.powf(d - 2.0);
        if val.abs() <= 1e-30 {
            0.0
        } else {
            val
        }
    }
}

fn lorenz_mie(
    x: f64,
    refrel: Complex<f64>,
    nstop: usize,
) -> Result<(Vec<Complex<f64>>, Vec<Complex<f64>>)> {
    let nmxx = 150_000;

    let y = refrel * x;
    let ymod = y.norm();
    let xstop = x + 4.0 * x.powf(1.0 / 3.0) + 2.0;
    let nmx = xstop.max(ymod) as usize + 15;
    if nmx > nmxx {
        bail!("ERROR: nmx > nmxx for |m|x = {}", ymod);
    }

    let mut d: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nmx];

    for n in 0..nmx - 1 {
        let en = nmx - 1 - n;
        let enr = (nmx - 1 - n) as f64;
        d[nmx - 1 - n] = (enr / y) - (1.0 / (d[en] + enr / y))
    }

    let mut psi_0 = x.cos();
    let mut psi_1 = x.sin();
    let mut chi_0 = -psi_1;
    let mut chi_1 = psi_0;
    let mut xi_1 = Complex::new(psi_1, -chi_1);

    let mut a: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nstop];
    let mut b: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nstop];

    for n in 0..nstop {
        let nr = n as f64;
        let psi = (2.0 * nr - 1.0) * psi_1 / x - psi_0;
        let chi = (2.0 * nr - 1.0) * chi_1 / x - chi_0;
        let xi = Complex::new(psi, -chi);
        a[n] = (d[n] / refrel + nr / x) * psi - psi_1;
        a[n] /= (d[n] / refrel + (nr / x)) * xi - xi_1;
        b[n] = (refrel * d[n] + nr / x) * psi - psi_1;
        b[n] /= (refrel * d[n] + nr / x) * xi - xi_1;
        psi_0 = psi_1;
        psi_1 = psi;
        chi_0 = chi_1;
        chi_1 = chi;
        xi_1 = Complex::new(psi_1, -chi_1);
    }
    Ok((a, b))
}

fn renormalize() {
    todo!()
}
