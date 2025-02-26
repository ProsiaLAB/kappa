//! Computes light scattering properties of randomly oriented
//! fractal dust aggregates by means of the modified mean field theory developed
//! in [Tazaki & Tanaka (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aac32d/meta). This code is also capable of computing the light
//! scattering solution based on the Rayleigh-Gans-Debye theory and the Mean field
//! theory.
//!
//! # References
//! Based on `OPFRACTAL v3.0` from `optool` by Ryo Tazaki.

use anyhow::bail;
use anyhow::Result;
use num_complex::Complex;

/// Complex vector type.
type ComplexVec = Vec<Complex<f64>>;

/// Computes light scattering properties of randomly oriented
/// fractal dust aggregates by means of the modified mean field theory developed
/// in [Tazaki & Tanaka (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aac32d/meta). This code is also capable of computing the light
/// scattering solution based on the Rayleigh-Gans-Debye theory and the Mean field
/// theory.
///
/// # Arguments
/// - `iqsca`:  Method switch for the light scatteing solver
///     - `iqsca = 1`: Rayleigh-Gans-Debye theory
///     - `iqsca = 2`: Mean field theory
///     - `iqsca = 3`: Modified mean-field theory
/// - `iqcor`: Switch for the two-point correction function
///     - `iqcor = 1`: Gaussian cut-off
///     - `iqcor = 2`: Exponential cut-off
///     - `iqcor = 3`: Fractal dimension cut-off
/// - `iqgeo`: Switch for the geometric cross-section
///     - `iqgeo = 1`: `π * rc^2`, where `rc` is the characteristic radius.
///     - `iqgeo = 2`: [Okuzumi et al. (2009)](https://iopscience.iop.org/article/10.1088/0004-637X/698/2/1122/meta)
///     - `iqgeo = 3`: [Tazaki (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.2811T/abstract)
/// - `nang`: Number of angles `(0, π/2)` for the scattering phase function.
/// - `pn`: Number of monomers in the aggregate.
/// - `r0`: Monomer radius.
/// - `df`: Fractal dimension.
/// - `k`: Prefactor for the fractal dimension.
/// - `lmd`: Wavelength of the incident light.
/// - `refrel`: Complex Refractive index of the monomer.
///
/// # Notes
/// - For `iqsca = 1`:
///   all outputs would be physically reasonable for the phase shift `< ~1`.
///
/// - For `iqsca = 2`:
///   The extinction cross section would be calculated without limitation.
///   However, the other outputs would be reliable for the phase shift < ~1.
///
/// - For `iqsca = 3`:
///   The extinction cross section would be calculated without limitation.
///   Scattering and absorption cross sections could be calculated
///   for the phase shift > 1, however too large phase shift may
///   cause some problem.
///   The asymmetry parameter and the sattering matrix elements would be
///   reliable for the phase shift < ~1.
///
/// - For the two-point correlation function, it is suggested to use `iqcor = 1` because
///   it is numerically stable and is likely to reproduce optical properties
///   of fractal aggregates.
///
/// # References
/// - `iqsca = 1`; `iqcor = 1`: [Tazaki et al. (2016), ApJ, 823, 70](https://iopscience.iop.org/article/10.3847/0004-637X/823/2/70/meta)
/// - `iqsca = 2`: [Botet et al. (1997), ApOpt, 36, 8791](https://opg.optica.org/abstract.cfm?uri=ao-36-33-8791)
/// - `iqsca = 3`: [Tazaki & Tanaka (2018), ApJ, 860, 79](https://iopscience.iop.org/article/10.3847/1538-4357/aac32d/meta)
/// - `iqcor = 2`: [Berry & Percival (1986), AcOpt, 33, 577](https://www.tandfonline.com/doi/abs/10.1080/713821987)
/// - `iqcor = 3`: [Botet et al. (1995), JPhA, 28, 297](https://iopscience.iop.org/article/10.1088/0305-4470/28/2/008)
/// - `iqgeo = 2`: [Okuzumi et al. (2009), ApJ, 707, 1247](https://iopscience.iop.org/article/10.1088/0004-637X/698/2/1122/meta)
/// - `iqgeo = 3`: [Tazaki (2021), MNRAS, 504, 2811](https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.2811T/abstract)
pub fn mean_scattering() {
    // maxwell_garnett_mixing();
    // structure_factor_integration();
    // lorenz_mie();
    todo!()
}

/// This function finds effective refractive index of an aggregates
/// by using Maxwell-Garnett mixing rule.
///
/// # Arguments
/// - `refrel`: Complex refractive index of the monomer.
/// - `f1`: Volume filling factor of the first material.
///
/// # Returns
/// - Effective refractive index of the aggregates.
///
/// # Notes
/// - eps_1 : Dielectric constant of the first material.
/// - eps_2 : Dielectric constant of vacuum.
/// - f1 : Volume filling factor of the first material.
fn maxwell_garnett_mixing(refrel: Complex<f64>, f1: f64) -> Complex<f64> {
    let eps_1 = refrel * refrel;
    let eps_2 = Complex::new(1.0, 0.0);
    let mg = eps_2 * (2.0 * f1 * (eps_1 - eps_2) + eps_1 + 2.0 * eps_2)
        / (eps_1 + 2.0 * eps_2 - f1 * (eps_1 - eps_2));

    mg.sqrt()
}

/// This subroutine performs integration of the static structure factor by
/// Botet et al. (1995) based on  `df`-dependent cut-off function for two-point
/// correlation function.
/// The statis structure factor `S(q)` is given by
///
///
///                     c*df     /x_max
///             S(q) =  ----  *  |  dx x^{df-2}*sin(q*Rg*x)*exp[-c*x^df],
///                     q*Rg     /0
///
///
/// where `x = u/Rg`; where `Rg` is the radius of gyration and `u` is the
/// distance between two monomers; `q` is the magnitude of scattering
/// vector, `df` is the fractal dimension, `c=0.5` is normalization factor.
///
/// Since we have a steep cut-off function for `x>1`, it would be enough
/// to take `x_max` slightly larger than unity.
/// Here we determine xmax so that `eta=c*x^df` becomes 25.
///
/// We use Simpson's rule for the numerical integration, and then,
/// the integration range `[0,xmax]` will be divided by `2*m` points.
///
/// # Arguments
/// - `d`: Fractal dimension.
/// - `q`: Magnitude of scattering vector.
/// - `r_g`: Radius of gyration.
///
/// # Returns
/// - Value of the static structure factor.
///
/// # Caution
/// For large `df(>2)` and large `qRg`, the numerical integration sometimes
/// returns a negative value of `S(q)`. Since `S(q)` is a power spectrum,
/// it must be positive. I also noticed that the same problem arises with
/// the code written by Pascal Rannou. Because of this ill behavior, I
/// prefer to use Gausian-cut-off model, which is more stable.
///
/// A not-too-bad tentative prescription is simply set `S(q)=0` when `S(q)<0`,
/// which means scattered waves from each monomer is incoherent. This is
/// what I did in if-sentence at the end of this subroutine.
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

/// Integrand of the static structure factor.
///
/// # Arguments
/// - `x`: Integration variable.
/// - `c`: Normalization factor.
/// - `d`: Fractal dimension.
/// - `q`: Magnitude of scattering vector.
/// - `r_g`: Radius of gyration.
///
/// # Returns
/// - Value of the integrand.
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

/// Calculate Lorenz-Mie scattering coefficients (an,bn) for a monomer particle.
///
/// Since monomer's size parameter is supposed not to be very large,
/// we use the simple Bohren & Huffman Mie algorithm.
///
/// # Arguments
/// - `x`: Size parameter.
/// - `refrel`: Complex refractive index of the monomer.
/// - `nstop`: Number of terms in the series expansion.
///
/// # Returns
/// - `a`: Scattering coefficient.
/// - `b`: Absorption coefficient.
///
/// # References
/// The original BHMIE code is taken from Bruce Draine's homepage:
/// <https://www.astro.princeton.edu/~draine/scattering.html>, with
/// slight modifications.
fn lorenz_mie(x: f64, refrel: Complex<f64>, nstop: usize) -> Result<(ComplexVec, ComplexVec)> {
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
