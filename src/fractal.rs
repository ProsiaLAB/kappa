//! Computes light scattering properties of randomly oriented
//! fractal dust aggregates by means of the modified mean field theory developed
//! in [Tazaki & Tanaka (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aac32d/meta). This code is also capable of computing the light
//! scattering solution based on the Rayleigh-Gans-Debye theory and the Mean field
//! theory.
//!
//! # References
//! Based on `OPFRACTAL v3.0` from `optool` by Ryo Tazaki.

use std::f64::consts::PI;

use anyhow::anyhow;
use anyhow::{Context, Result};
use extensions::types::{CMatrix, CTensor, CVector, RMatrix, RVector};
use ndarray::s;
use num_complex::{Complex, Complex64};

use crate::geofractal::get_geometric_cross_section_tazaki;
use crate::geofractal::{AFactor, IntegrationMethod};
use crate::utils::gamma::gamma;
use crate::utils::{bessel, legendre, linalg, special};

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
    pub refrel: Complex64,
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
    pub ang: RVector,
    pub smat: RMatrix,
    pub pf: RVector,
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
pub fn mean_scattering(fracc: &FractalConfig) -> Result<FractalResult> {
    let k = 2.0 * PI / fracc.lmd; // Wavenumber
    let r_g = fracc.r0 * (fracc.pn / fracc.k0).powf(1.0 / fracc.df); // Radius of gyration of the aggregate
    let r_c = (5.0f64 / 3.0f64).powf(r_g); // Characteristic radius of the aggregate
    let x_g = k * r_g; // Size parameter of the aggregate
    let x0 = k * fracc.r0; // Size parameter of the monomer

    let eta = 25.0;
    let q_r_g_critical = 26.0;

    // Truncation order of the monomer's scattering field
    let xmax = x0 + 4.0 * x0.powf(1.0 / 3.0) + 2.0;
    let nmax = xmax as usize;

    if fracc.nang <= 1 {
        return Err(anyhow!("LowScatteringAngleResolution"));
    }

    if fracc.pn < 1.0 {
        return Err(anyhow!("NotEnoughMonomers"));
    }

    if fracc.df > 3.0 {
        return Err(anyhow!("ExceedsMaxFractalDimension"));
    }

    if (2 * nmax) >= 500 {
        eprintln!("WARNING: The truncation order of monomer's scattered light");
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

    let (an, bn) = lorenz_mie(x0, fracc.refrel, nmax)?;

    let mut ad = CMatrix::zeros((nmax, 2));
    let mut dd = CMatrix::zeros((nmax, 2));

    ad.slice_mut(s![.., 0]).assign(&an);
    ad.slice_mut(s![.., 1]).assign(&bn);

    // Solve multiple scattering?
    match fracc.solver {
        FractalSolver::RayleighGansDebye => {
            dd.slice_mut(s![.., 0]).assign(&ad.slice(s![.., 0]));
            dd.slice_mut(s![.., 1]).assign(&ad.slice(s![.., 1]));
        }
        FractalSolver::MeanField | FractalSolver::ModifiedMeanField => {
            let r = mean_field(fracc, &ad, x_g, nmax)?;
            for n in 0..nmax {
                dd[[n, 0]] = r[2 * n - 1];
                dd[[n, 1]] = r[2 * n];
            }
        }
    }

    let dang = PI / 2.0 / (fracc.nang - 1) as f64;
    let (s1, s2) = renormalize(fracc, &dd, dang, nmax);

    let mut smat = RMatrix::zeros((2 * fracc.nang, 4));
    let ang = (0..(2 * fracc.nang))
        .map(|j| j as f64 * dang)
        .collect::<RVector>();

    for j in 0..(2 * fracc.nang) {
        let q = 2.0 * k * (0.5 * ang[j]).sin();
        let s11 = 0.5 * (s2[j].norm().powi(2) + s1[j].norm().powi(2));
        let s12 = 0.5 * (s2[j].norm().powi(2) - s1[j].norm().powi(2));
        let s33 = (s1[j] * s2[j].conj()).re;
        let s34 = (s2[j] * s1[j].conj()).im;
        let sq = match fracc.cutoff {
            FractalCutoff::Gaussian => {
                if fracc.df == 3.0 {
                    if q * q * r_g * r_g >= eta {
                        0.0
                    } else {
                        (-q * q * r_g * r_g / 3.0).exp()
                    }
                } else {
                    let al = 0.5 * fracc.df;
                    let bb = 1.5;
                    let xx = -q * q * r_g * r_g / fracc.df;
                    if q * r_g < q_r_g_critical {
                        special::confluent_hypergeometric(al, bb, xx)
                    } else {
                        let cc = PI.sqrt() * fracc.df.powf(0.5 * fracc.df)
                            / (2.0 * gamma(0.5 * (3.0 - fracc.df)));
                        cc * (q * r_g).powf(-fracc.df)
                    }
                }
            }
            FractalCutoff::Exponential => {
                let xi = (2.0 / (fracc.df * (fracc.df + 1.0))).sqrt() * r_g;
                if j == 0 {
                    1.0
                } else {
                    ((fracc.df - 1.0) * (q * xi).atan()).sin()
                        / (q * xi
                            * (fracc.df - 1.0)
                            * (1.0 + (q * xi).powi(2)).powf((fracc.df - 1.0) / 2.0))
                }
            }
            FractalCutoff::FractalDimension => structure_factor_integration(fracc.df, q, r_g),
        };

        smat[[j, 0]] = fracc.pn * s11 * (1.0 + (fracc.pn - 1.0) * sq);
        smat[[j, 1]] = fracc.pn * s12 * (1.0 + (fracc.pn - 1.0) * sq);
        smat[[j, 2]] = fracc.pn * s33 * (1.0 + (fracc.pn - 1.0) * sq);
        smat[[j, 3]] = fracc.pn * s34 * (1.0 + (fracc.pn - 1.0) * sq);

        if smat[[j, 0]] < 0.0 {
            return Err(anyhow!("IllegalMatrixElement"));
        }
    }

    let mut c_sca: f64 = (0..(fracc.nang - 1))
        .map(|j| {
            dang * (smat[[j, 0]] * (ang[j].sin())
                + 4.0 * smat[[j + 1, 0]] * (ang[j + 1].sin())
                + smat[[j + 2, 0]] * (ang[j + 2].sin()))
                / 3.0
        })
        .sum::<f64>()
        * 2.0
        * PI
        / (k * k);

    let pf = smat.slice(s![.., 0]).to_owned() / c_sca / (k * k);
    let g = (0..(fracc.nang - 1))
        .map(|j| {
            dang * (pf[j] * (ang[j].sin()) * (ang[j].cos())
                + 4.0 * pf[j + 1] * (ang[j + 1].sin()) * (ang[j + 1].cos())
                + pf[j + 2] * (ang[j + 2].sin()) * (ang[j + 2].cos()))
                / 3.0
        })
        .sum::<f64>()
        * 2.0
        * PI;
    let nrm = (0..(fracc.nang - 1))
        .map(|j| {
            dang * (pf[j] * (ang[j].sin())
                + 4.0 * pf[j + 1] * (ang[j + 1].sin())
                + pf[j + 2] * (ang[j + 2].sin()))
                / 3.0
        })
        .sum::<f64>()
        * 2.0
        * PI;

    // Check normalization of phase function
    if (nrm - 1.0).abs() > 1e-3 {
        return Err(anyhow!("UnnormalizedPhaseFunction"));
    }
    let c_abs: f64;
    let c_ext: f64;
    (c_sca, c_abs, c_ext) = match fracc.solver {
        FractalSolver::RayleighGansDebye => {
            let c_abs = (0..nmax)
                .map(|j| {
                    let jf = j as f64;
                    let cn1 = 1.0 / ad[[j, 0]].conj() - 1.0;
                    let cn2 = 1.0 / ad[[j, 1]].conj() - 1.0;

                    (2.0 * jf + 3.0)
                        * (cn1 * ad[[j, 0]].norm() * ad[[j, 0]].norm()
                            + cn2 * ad[[j, 1]].norm() * ad[[j, 1]].norm())
                        .re
                })
                .sum::<f64>()
                * 2.0
                * PI
                * fracc.pn
                / k
                / k;
            let c_ext = c_abs + c_sca;
            (c_sca, c_abs, c_ext)
        }
        FractalSolver::MeanField => {
            let c_ext = (0..nmax)
                .map(|j| {
                    let jf = j as f64;
                    (2.0 * jf + 3.0) * (dd[[j, 0]] + dd[[j, 1]]).re
                })
                .sum::<f64>()
                * 2.0
                * PI
                * fracc.pn
                / k
                / k;
            let c_abs = c_ext - c_sca;
            (c_sca, c_abs, c_ext)
        }
        FractalSolver::ModifiedMeanField => {
            let c_ext = (0..nmax)
                .map(|j| {
                    let jf = j as f64;
                    (2.0 * jf + 3.0) * (dd[[j, 0]] + dd[[j, 1]]).re
                })
                .sum::<f64>()
                * 2.0
                * PI
                * fracc.pn
                / k
                / k;
            let mut c_abs = (0..nmax)
                .map(|j| {
                    let cn1 = 1.0 / ad[[j, 0]].conj() - 1.0;
                    let cn2 = 1.0 / ad[[j, 1]].conj() - 1.0;
                    (2.0 * j as f64 + 3.0)
                        * (cn1 * ad[[j, 0]].norm() * ad[[j, 0]].norm()
                            + cn2 * ad[[j, 1]].norm() * ad[[j, 1]].norm())
                        .re
                })
                .sum::<f64>()
                * 2.0
                * PI
                * fracc.pn
                / k
                / k;

            let gc = match fracc.geometry {
                FractalGeometry::Circular => PI * r_c * r_c,
                FractalGeometry::Okuzumi => {
                    get_geometric_cross_section_okuzumi(fracc.pn, fracc.r0, r_c)
                }
                FractalGeometry::Tazaki => {
                    get_geometric_cross_section_tazaki(
                        &IntegrationMethod::AnAn,
                        &AFactor::Unity,
                        &fracc.cutoff,
                        fracc.pn,
                        fracc.k0,
                        fracc.df,
                    )? * fracc.pn
                        * PI
                        * fracc.r0
                        * fracc.r0
                }
            };
            let tau = c_abs / gc;
            if tau >= 10.0 {
                c_abs = gc;
            } else {
                c_abs = gc * (1.0 - (-tau).exp());
            }
            c_abs = c_abs.max(c_ext - c_sca);
            c_sca = c_ext - c_abs;
            (c_sca, c_abs, c_ext)
        }
    };

    let fracr = FractalResult {
        c_ext,
        c_sca,
        c_abs,
        g,
        dphi,
        ang,
        smat,
        pf,
    };
    Ok(fracr)
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
/// - `eps_1` : Dielectric constant of the first material.
/// - `eps_2` : Dielectric constant of vacuum.
/// - `f1` : Volume filling factor of the first material.
fn maxwell_garnett_mixing(refrel: Complex64, f1: f64) -> Complex64 {
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
/// $$
/// S(q) = \frac{c \cdot df}{q R_g} \int_0^{x_{\max}} dx \, x^{df-2} \sin(q R_g x) e^{-c x^{df}}
/// $$
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
fn lorenz_mie(x: f64, refrel: Complex64, nstop: usize) -> Result<(CVector, CVector)> {
    let nmxx = 150_000;

    let y = refrel * x;
    let ymod = y.norm();
    let xstop = x + 4.0 * x.powf(1.0 / 3.0) + 2.0;
    let nmx = xstop.max(ymod) as usize + 15;
    if nmx > nmxx {
        return Err(anyhow!("MieOverflow"));
    }

    let mut d = CVector::zeros(nmx);

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

    let mut a = CVector::zeros(nstop);
    let mut b = CVector::zeros(nstop);

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
fn renormalize(fracc: &FractalConfig, d: &CMatrix, dang: f64, nmax: usize) -> (CVector, CVector) {
    let amu = (0..fracc.nang)
        .map(|j| j as f64 * dang)
        .map(|theta| theta.cos())
        .collect::<RVector>();

    let mut pi = RVector::zeros(fracc.nang);
    let mut pi0 = RVector::zeros(fracc.nang);
    let mut pi1 = RVector::zeros(fracc.nang);

    let mut tau = RVector::zeros(fracc.nang);

    let nn = 2 * fracc.nang;
    let mut s1 = CVector::zeros(nn);
    let mut s2 = CVector::zeros(nn);

    let mut p = -1.0;
    // let mut an = Complex::new(0.0, 0.0);
    // let mut an1 = Complex::new(0.0, 0.0);
    // let mut bn = Complex::new(0.0, 0.0);
    // let mut bn1 = Complex::new(0.0, 0.0);

    for n in 0..nmax {
        let nf = n as f64;
        let f_n = (2.0 * nf + 3.0) / (nf + 1.0) / (nf + 2.0);
        // if n > 1 {
        //     an1 = an;
        //     bn1 = bn;
        // }
        let an = d[[n, 0]];
        let bn = d[[n, 1]];

        for j in 0..fracc.nang {
            pi[j] = pi1[j];
            tau[j] = nf * amu[j] * pi[j] - (nf + 2.0) * pi0[j];
            s1[j] += f_n * (an * pi[j] + bn * tau[j]);
            s2[j] += f_n * (an * tau[j] + bn * pi[j]);
        }
        p = -p;
        for j in 0..(fracc.nang - 1) {
            let jj = 2 * fracc.nang - j - 1;
            s1[jj] += f_n * p * (an * pi[j] - bn * tau[j]);
            s2[jj] += f_n * p * (-an * tau[j] + bn * pi[j]);
        }

        for j in 0..fracc.nang {
            pi1[j] = ((2.0 * nf + 3.0) * amu[j] * pi[j] - (nf + 2.0) * pi0[j]) / (nf + 1.0);
            pi0[j] = pi[j];
        }
    }
    (s1, s2)
}

/// Solves the mean field theory.
///
/// # PART I.
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
///
/// # PART II.
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
fn mean_field(fracc: &FractalConfig, ad: &CMatrix, x_g: f64, nmax: usize) -> Result<CVector> {
    let jm = 400; // Number of grid points for Gauss-Legendre quadrature

    // Integration limits
    let x1 = -1.0;
    let x2 = 1.0;

    // Gauss-Legendre quadrature
    let (x, w) = legendre::gauss_legendre(x1, x2, jm);

    let order = 1;
    let pmax = 2 * nmax;

    // Store values of Legendre polynomials (P_n^m) and their derivatives
    let mut al1n = RMatrix::zeros((jm, nmax));
    let mut ln = RMatrix::zeros((jm, pmax));
    let mut dln = RMatrix::zeros((jm, pmax));

    for (j, xj) in x.iter().enumerate() {
        let (pmn, _) =
            legendre::lpmns(order, nmax, *xj).context("IntegrationFailure: In Legendre")?;
        let (lp, dlp) = legendre::lpn(pmax, *xj);
        al1n.slice_mut(s![j, ..]).assign(&pmn);
        ln.slice_mut(s![j, ..]).assign(&lp);
        dln.slice_mut(s![j, ..]).assign(&dlp);
    }

    // Pre-compute the spherical Bessel functions and the T-matrix
    let mut s_p = CVector::zeros(pmax);

    for (p, val) in s_p.iter_mut().enumerate() {
        *val = bessel::int_sph_bessel(fracc, x_g, p).context("IntegrationFailure: In Bessel")?;
    }

    let mut t = CTensor::zeros((2, nmax, nmax));
    for nu in 0..nmax {
        let nuf = nu as f64;
        for n in 0..nmax {
            let nf = n as f64;
            let pmin = nu.abs_diff(n);
            let pmax = n + nu;
            let mut sum_a = Complex::new(0.0, 0.0);
            let mut sum_b = Complex::new(0.0, 0.0);
            for (p, val) in s_p.iter().enumerate().take(pmax + 1).skip(pmin) {
                let val = *val;
                let pf = p as f64;
                let mut anunp = 0.0;
                let mut bnunp = 0.0;
                if pmax % 2 == p % 2 {
                    anunp = w
                        .iter()
                        .zip(al1n.row(nu))
                        .zip(al1n.row(n))
                        .zip(ln.row(p))
                        .map(|(((w, al1n_nu), al1n_n), ln_p)| w * al1n_nu * al1n_n * ln_p)
                        .sum();
                } else {
                    bnunp = w
                        .iter()
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

    let mut s = CMatrix::ones((pmax, 2));
    let mut y = CVector::zeros(pmax);

    for n in 0..nmax {
        for nu in 0..nmax {
            s[[2 * n - 1, 2 * n - 1]] = (fracc.pn - 1.0) * ad[[0, n]] * t[[0, nu, n]];
            s[[2 * n - 1, 2 * n]] = (fracc.pn - 1.0) * ad[[0, n]] * t[[1, nu, n]];
            s[[2 * n, 2 * n - 1]] = (fracc.pn - 1.0) * ad[[1, n]] * t[[0, nu, n]];
            s[[2 * n, 2 * n]] = (fracc.pn - 1.0) * ad[[1, n]] * t[[1, nu, n]];
        }
        y[2 * n - 1] = ad[[0, n]];
        y[2 * n] = ad[[1, n]];
    }

    // Invert the T-matrix elements S
    // and find mean-field vector r
    let r = linalg::complex_matrix_inverse(pmax, 2 * pmax, &y, &s).context("SingularMatrix")?;

    Ok(r)
}

fn get_geometric_cross_section_okuzumi(pn: f64, a0: f64, ac: f64) -> f64 {
    let k0_bcca = 1.04; // Tazaki+2016; BCCA with oblique collisions
    let df_bcca = 1.90; // Tazaki+2016; BCCA with oblique collisions

    let ac_bcca = (5.0f64 / 3.0f64).sqrt() * a0 * (pn / k0_bcca).powf(1.0 / df_bcca);
    let gc0 = PI * a0 * a0;
    let gcc = PI * ac * ac;

    let gc_bcca2 = PI * ac_bcca * ac_bcca;

    let gc_bcca = if pn <= 16.0 {
        12.5 * pn.powf(0.685) * (-2.53 / pn.powf(0.0920)).exp() * gc0
    } else {
        0.352 * pn + 0.566 * pn.powf(0.862) * gc0
    };

    (1.0 / (1.0 / gc_bcca + 1.0 / gcc - 1.0 / gc_bcca2)).min(pn * PI * a0 * a0)
}
