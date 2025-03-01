//! Computes light scattering properties of randomly oriented
//! fractal dust aggregates by means of the modified mean field theory developed
//! in [Tazaki & Tanaka (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aac32d/meta). This code is also capable of computing the light
//! scattering solution based on the Rayleigh-Gans-Debye theory and the Mean field
//! theory.
//!
//! # References
//! Based on `OPFRACTAL v3.0` from `optool` by Ryo Tazaki.

use std::f64::consts::PI;

use anyhow::Result;
use ndarray::s;
use ndarray::Array;
use ndarray::Array2;
use ndarray::Array3;
use num_complex::Complex;

use crate::utils::bessel;
use crate::utils::legendre;

/// Complex vector type.
type ComplexVec = Vec<Complex<f64>>;

pub struct FractalConfig {
    pub solver: FractalSolver,
    pub cutoff: FractalCutoff,
    pub geometry: FractalGeometry,
    pub nang: usize,
    pub pn: f64,
    pub r0: f64,
    pub df: f64,
    pub k0: f64,
    pub lmd: f64,
    pub refrel: Complex<f64>,
    pub gauss_leg_grid_size: usize,
    pub truncation_order: usize,
}

impl Default for FractalConfig {
    fn default() -> Self {
        FractalConfig {
            solver: FractalSolver::ModifiedMeanField,
            cutoff: FractalCutoff::Gaussian,
            geometry: FractalGeometry::Circular,
            nang: 100,
            pn: 1.0,
            r0: 0.1,
            df: 2.0,
            k0: 1.0,
            lmd: 0.55,
            refrel: Complex::new(1.5, 0.0),
            gauss_leg_grid_size: 400,
            truncation_order: 500,
        }
    }
}

impl FractalConfig {
    pub fn new() -> Self {
        Self::default()
    }
}

pub struct FractalResult {
    pub c_ext: f64,
    pub c_sca: f64,
    pub c_abs: f64,
    pub g: f64,
    pub dphi: f64,
    pub ang: Vec<f64>,
    pub smat: Array2<f64>,
    pub pf: Vec<f64>,
}

pub enum FractalError {
    ScatteringAngleResolutionTooLow,
    InsufficientNumberOfMonomers,
    FractalDimensionTooLarge,
    MieOverflow,
    IntegrationFailure(String),
}

pub enum FractalSolver {
    RayleighGansDebye,
    MeanField,
    ModifiedMeanField,
}

pub enum FractalCutoff {
    Gaussian,
    Exponential,
    FractalDimension,
}

#[derive(PartialEq)]
pub enum FractalGeometry {
    Circular,
    Okuzumi,
    Tazaki,
}

pub struct TMatrixConfig {
    pub numax: usize,
    pub nmax: usize,
    pub jm: usize,
    pub al1n: Array2<f64>,
    pub ln: Array2<f64>,
    pub dln: Array2<f64>,
    pub w: Vec<f64>,
    pub sp: Vec<Complex<f64>>,
}

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
pub fn mean_scattering(fracc: &FractalConfig) -> Result<(), FractalError> {
    let jm = 400; // Number of grid points for Gauss-Legendre quadrature

    let k = 2.0 * PI / fracc.lmd; // Wavenumber
    let r_g = fracc.r0 * (fracc.pn / fracc.k0).powf(1.0 / fracc.df); // Radius of gyration of the aggregate
    let r_c = (5.0f64 / 3.0f64).powf(r_g); // Characteristic radius of the aggregate
    let x_g = k * r_g; // Size parameter of the aggregate
    let x0 = k * fracc.r0; // Size parameter of the monomer

    // Truncation order of the monomer's scattering field
    let x_stop = x0 + 4.0 * x0.powf(1.0 / 3.0) + 2.0;
    let nstop = x_stop as usize;
    let numax = nstop;
    let nmax = nstop;

    if fracc.nang <= 1 {
        return Err(FractalError::ScatteringAngleResolutionTooLow);
    }

    if fracc.pn < 1.0 {
        return Err(FractalError::InsufficientNumberOfMonomers);
    }

    if fracc.df > 3.0 {
        return Err(FractalError::FractalDimensionTooLarge);
    }

    if numax + nmax >= 500 {
        eprintln!("WARNING: numax + nmax >= 500");
        eprintln!("         The truncation order of monomer's scattered light");
        eprintln!("         field exceeds the maximum value (=500).          ");
        eprintln!("         This may cause a code crush at computations of   ");
        eprintln!("         the spherical Bessel function of the first kind. ");
    }

    // Compute the phase shift (cf: Eq. 9 in Tazaki & Tanaka 2018)
    let ff = fracc.pn * (fracc.r0 / r_c).powf(3.0);
    let mgmref = maxwell_garnett_mixing(fracc.refrel, ff);
    let dphic = 2.0 * k * r_c * (mgmref - 1.0).norm();
    let dphi0 = 2.0 * x0 * (fracc.refrel - 1.0).norm();
    let dphi = dphic.max(dphi0);

    if dphi >= 1.0 {
        eprintln!("WARNING: dphi >= 1");
        eprintln!("         The phase shift by an aggregate exceeds unity.");
        eprintln!("         Output of scattering matrix elements are not  ");
        eprintln!("         physically guaranteed.                        ");
    }

    let (an, bn) = lorenz_mie(x0, fracc.refrel, nstop)?;

    let mut ad: Array2<Complex<f64>> = Array2::zeros((nstop, 2));
    let mut dd: Array2<Complex<f64>> = Array2::zeros((nstop, 2));

    ad.slice_mut(s![.., 0]).assign(&Array::from_vec(an));
    ad.slice_mut(s![.., 1]).assign(&Array::from_vec(bn));

    // Solve multiple scattering?
    match fracc.solver {
        FractalSolver::RayleighGansDebye => {
            dd.slice_mut(s![.., 0]).assign(&ad.slice(s![.., 0]));
            dd.slice_mut(s![.., 1]).assign(&ad.slice(s![.., 1]));
        }
        FractalSolver::MeanField | FractalSolver::ModifiedMeanField => {
            mean_field(fracc, x_g, jm, nstop, nmax, numax)?;
        }
    }

    // let fracr: FractalResult;
    Ok(())
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
            let ir = i as f64;
            let x0 = 2.0 * ir * h;
            let x1 = (2.0 * ir + 1.0) * h;
            let x2 = (2.0 * ir + 2.0) * h;
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
fn lorenz_mie(
    x: f64,
    refrel: Complex<f64>,
    nstop: usize,
) -> Result<(ComplexVec, ComplexVec), FractalError> {
    let nmxx = 150_000;

    let y = refrel * x;
    let ymod = y.norm();
    let xstop = x + 4.0 * x.powf(1.0 / 3.0) + 2.0;
    let nmx = xstop.max(ymod) as usize + 15;
    if nmx > nmxx {
        return Err(FractalError::MieOverflow);
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
        let nf = n as f64;
        let psi = (2.0 * nf - 1.0) * psi_1 / x - psi_0;
        let chi = (2.0 * nf - 1.0) * chi_1 / x - chi_0;
        let xi = Complex::new(psi, -chi);
        a[n] = (d[n] / refrel + nf / x) * psi - psi_1;
        a[n] /= (d[n] / refrel + (nf / x)) * xi - xi_1;
        b[n] = (refrel * d[n] + nf / x) * psi - psi_1;
        b[n] /= (refrel * d[n] + nf / x) * xi - xi_1;
        psi_0 = psi_1;
        psi_1 = psi;
        chi_0 = chi_1;
        chi_1 = chi;
        xi_1 = Complex::new(psi_1, -chi_1);
    }
    Ok((a, b))
}

/// Calculate scattering amplitude S1, S2 from given scattering coefficients d.
///
/// # References
/// The original BHMIE code is taken from Bruce Draine's homepage:
/// <https://www.astro.princeton.edu/~draine/scattering.html>, with
/// slight modifications.
fn renormalize() {}

/// Solves the mean field theory.
///
/// Calculate a(nu,n,p) and b(nu,n,p):
///
///                   2p + 1  /+1
///       a(nu,n,p) = -------- | dx P_nu^1(x) * P_n^1(x) * P_p(x),
///                      2     /-1
///     
///                    2p + 1  /+1                         dP_p(x)
///       b(nu,n,p) = -------- | dx P_nu^1(x) * P_n^1(x) * -------,
///                      2     /-1                            dx
///
/// where `P_n^m` is the associated Legendre function (`n`: degree, `m`: order),
///  `P_n` is the Legendre polynominal function
///  (see Equations (29, 30) in Tazaki & Tanaka 2018).
///  The integration is performed with the Gauss-Legendre quadrature.
fn mean_field(
    fracc: &FractalConfig,
    x_g: f64,
    jm: usize,
    nstop: usize,
    nmax: usize,
    numax: usize,
) -> Result<(), FractalError> {
    // Solve the mean field theory
    let x1 = -1.0;
    let x2 = 1.0;

    let (x, w) = legendre::gauss_legendre(x1, x2, jm);

    let order = 1;
    let degmax = nstop;
    let pmax = 2 * nstop;

    let mut al1n: Array2<f64> = Array2::zeros((jm, degmax));
    let mut ln: Array2<f64> = Array2::zeros((jm, pmax));
    let mut dln: Array2<f64> = Array2::zeros((jm, pmax));

    for (j, xj) in x.iter().enumerate().take(jm) {
        let (pmn, _) = legendre::lpmns(order, degmax, *xj)
            .map_err(|_| FractalError::IntegrationFailure("In Legendre".to_string()))?;
        let (lp, dlp) = legendre::lpn(pmax, *xj);
        al1n.slice_mut(s![j, ..]).assign(&Array::from_vec(pmn));
        ln.slice_mut(s![j, ..]).assign(&Array::from_vec(lp));
        dln.slice_mut(s![j, ..]).assign(&Array::from_vec(dlp));
    }

    let mut s_p = vec![Complex::new(0.0, 0.0); numax + nmax];

    for (p, val) in s_p.iter_mut().enumerate().take(numax + nmax) {
        *val = bessel::int_sph_bessel(fracc, x_g, p).map_err(|_| {
            FractalError::IntegrationFailure(
                "Something went wrong in Bessel integration".to_string(),
            )
        })?;
    }

    get_translation_matrix_coefficients(numax, nmax, jm, &al1n, &ln, &dln, &w, &s_p);

    Ok(())
}

/// Calculate translation matrix coefficients
///
/// Translation coefficients: A and B (Equation # from Tazaki & Tanaka 2018)
/// T(1,nu,n) : \bar{A}_{1,n}^{1,nu} defined in Equation (14)
/// T(2,nu,n) : \bar{A}_{1,n}^{1,nu} defined in Equation (15)
/// anunp     : a(nu,n,p)            defined in Equation (29)
/// bnunp     : b(nu,n,p)            defined in Equation (30)
///
/// Note about the rule of rum with respect to the index p
/// (RT and Botet, in private communication, 2017).
///
/// Just before Equation (4) in Botet et al. (1997), it is written that
/// p should have the same parity as n+nu. However, this is only true
/// true for a(nu,n,p), and not for b(nu,n,p)
/// Thus, in this code, index p runs ALL integers between |nu-n| and nu+n.
///
/// Tips for efficient computation.
/// The parity of the integrand of a and b leads following properties:
/// - a(nu,n,p) .ne. 0 and b(nu,n,p) = 0 when p has the same parity as n+nu
/// - b(nu,n,p) .ne. 0 and a(nu,n,p) = 0 when p has the opposite parity to n+nu
///
///--------------------------------------------------------------------------------
fn get_translation_matrix_coefficients(
    numax: usize,
    nmax: usize,
    jm: usize,
    al1n: &Array2<f64>,
    ln: &Array2<f64>,
    dln: &Array2<f64>,
    w: &[f64],
    sp: &[Complex<f64>],
) -> Array3<Complex<f64>> {
    let mut t: Array3<Complex<f64>> = Array3::zeros((2, numax, nmax));
    for nu in 0..numax {
        let nuf = nu as f64;
        for n in 0..nmax {
            let nf = n as f64;
            let pmin = nu.abs_diff(n);
            let pmax = n + nu;
            let mut sum_a = Complex::new(0.0, 0.0);
            let mut sum_b = Complex::new(0.0, 0.0);
            for (p, val) in sp.iter().enumerate().take(pmax + 1).skip(pmin) {
                let val = *val;
                let pf = p as f64;
                let mut anunp = 0.0;
                let mut bnunp = 0.0;
                if pmax % 2 == p % 2 {
                    anunp = w
                        .iter()
                        .take(jm)
                        .zip(al1n.row(nu))
                        .zip(al1n.row(n))
                        .zip(ln.row(p))
                        .map(|(((w, al1n_nu), al1n_n), ln_p)| w * al1n_nu * al1n_n * ln_p)
                        .sum();
                } else {
                    bnunp = w
                        .iter()
                        .take(jm)
                        .zip(al1n.row(nu))
                        .zip(al1n.row(n))
                        .zip(dln.row(p))
                        .map(|(((w, al1n_nu), al1n_n), dln_p)| w * al1n_nu * al1n_n * dln_p)
                        .sum();
                }
                anunp *= 0.5 * (2.0 * pf + 1.0);
                bnunp *= 0.5 * (2.0 * pf + 1.0);

                sum_a += (nf * (nf + 1.0) + nuf * (nuf + 1.0) - pf * (pf + 1.0)) * anunp * val;
                sum_b += bnunp * val;
            }
            t[[0, nu, n]] = sum_a * (2.0 * nuf + 1.0) / (nf * (nf + 1.0) * nuf * (nuf + 1.0));
            t[[1, nu, n]] = sum_b * 2.0 * (2.0 * nuf + 1.0) / (nf * (nf + 1.0) * nuf * (nuf + 1.0));
        }
    }
    t
}
