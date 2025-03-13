//! The heart of `kappa`.
//!
//!

use std::f64::consts::PI;
use std::mem::swap;

use anyhow::Result;
use num_complex::ComplexFloat;
use num_complex::{Complex, Complex64};

pub struct KappaConfig {
    pub amin: f64,
    pub amax: f64,
    pub apow: f64,
    pub amean: f64,
    pub asigma: f64,
    pub na: usize,
    pub sizedis: SizeDistribution,
    pub lmin: f64,
    pub lmax: f64,
    pub nlam: usize,
    pub nang: usize,
    pub chop_angle: f64,
    pub pcore: f64,
    pub pmantle: f64,
    pub nmat: usize,
    pub nmant: usize,
    pub method: KappaMethod,
    pub fmax: f64,
    pub write_fits: bool,
    pub write_scatter: bool,
    pub for_radmc: bool,
    pub nsparse: usize,
    pub mmf_struct: f64,
    pub mmf_a0: f64,
    pub split: bool,
    pub blend_only: bool,
}

impl Default for KappaConfig {
    fn default() -> Self {
        KappaConfig {
            amin: 0.05,
            amax: 3000.0,
            apow: 3.5,
            amean: 0.0,
            asigma: 0.0,
            na: 0,
            sizedis: SizeDistribution::Apow,
            lmin: 0.05,
            lmax: 1e4,
            nlam: 300,
            nang: 180,
            chop_angle: 0.0,
            pcore: 0.0,
            pmantle: 0.0,
            nmat: 0,
            nmant: 0,
            method: KappaMethod::DHS,
            fmax: 0.8,
            write_fits: false,
            write_scatter: false,
            for_radmc: false,
            nsparse: 0,
            mmf_struct: 0.0,
            mmf_a0: 0.0,
            split: true,
            blend_only: false,
        }
    }
}

pub trait SpecialConfigs {
    fn diana() -> Self;
    fn dsharp() -> Self;
    fn dsharp_no_ice() -> Self;
}

impl SpecialConfigs for KappaConfig {
    fn diana() -> Self {
        KappaConfig::default()
    }

    fn dsharp() -> Self {
        KappaConfig::default()
    }

    fn dsharp_no_ice() -> Self {
        KappaConfig::default()
    }
}

#[derive(PartialEq)]
pub enum SizeDistribution {
    Apow,
    File,
    Norm,
    Lognorm,
}

impl SizeDistribution {
    fn validate(&self, kpc: &KappaConfig) -> Result<(), KappaError> {
        match self {
            SizeDistribution::Apow if kpc.apow < 0.0 => Err(KappaError::UnexpectedApow),
            SizeDistribution::Norm => {
                if kpc.amean <= 0.0 {
                    return Err(KappaError::InvalidSizeParam(
                        "`amean` must be positive for (log-)normal distribution".to_string(),
                    ));
                }
                if kpc.asigma == 0.0 {
                    return Err(KappaError::InvalidSizeParam(
                        "`asigma` cannot be zero for (log-)normal distribution".to_string(),
                    ));
                }
                if kpc.asigma > 0.0 {
                    return Err(KappaError::ForceLogNormal);
                }
                Ok(())
            }
            _ => Ok(()),
        }
    }
}

pub enum KappaMethod {
    DHS,
    MMF,
    CDE,
}

#[derive(Debug)]
pub enum KappaError {
    TooManyMaterials,
    NotEnoughMaterials,
    ZeroDensity,
    MissingDensity,
    InvalidMaterialKey,
    InvalidLNKFile,
    DuplicateDensity,
    MaterialBeforeStandard,
    InvalidSizeInput,
    InvalidSizeParam(String),
    SizeDistFileNotFound,
    DeltaExceedsSize,
    InvalidArgument,
    InvalidWavelengthInput,
    MissingSizeParamLimit,
    MissingSizeParamDHSLimit,
    TooManySparseFileRanges,
    NotEnoughWavelengthsForSparse,
    MissingFeatureValue,
    InvalidFeature,
    InvalidPorosity,
    IncompatibleArguments(String),
    UnnormalizedAbundance,
    SetSampling,
    SamplingRequired,
    SizeReorderRange,
    SizeForceUnitSampling,
    WavelengthReorderRange,
    WavelengthForceUnitSampling,
    UnexpectedApow,
    ForceLogNormal,
    DisableSplit,
}

/// Mueller matrix structure
pub struct Mueller {
    pub f11: Vec<f64>,
    pub f12: Vec<f64>,
    pub f22: Vec<f64>,
    pub f33: Vec<f64>,
    pub f44: Vec<f64>,
    pub f34: Vec<f64>,
}

/// Particle
pub struct Particle {
    pub rv: f64,
    pub rvmin: f64,
    pub rvmax: f64,
    pub rho: f64,
    pub k_abs: Vec<f64>,
    pub k_sca: Vec<f64>,
    pub k_ext: Vec<f64>,
    pub g: Vec<f64>,
    pub f: Vec<Mueller>,
    pub trust: Vec<bool>,
    pub is_ok: bool,
    pub is_ok_lmin: f64,
}

/// Defines a material component which is initialized
/// at runtime by consuming a LNK file.
///
/// See [`crate::components::StaticComponent`] for an equivalent
/// version which can is instantiated statically at compile time.
#[derive(Debug)]
pub struct Component {
    /// Name of the component.
    pub name: String,
    /// Class of the component.
    pub class: String,
    /// State of the component.
    pub state: String,
    /// Density of the component.
    pub rho: f64,
    /// size
    pub size: usize,
    /// Wavelengths.
    pub l0: Vec<f64>,
    /// Refractive index.
    pub n0: Vec<f64>,
    /// Extinction coefficient.
    pub k0: Vec<f64>,
}

pub fn run(kpc: &mut KappaConfig) -> Result<(), KappaError> {
    prepare_inputs(kpc)?;

    compute_kappa(kpc)?;
    Ok(())
}

fn prepare_inputs(kpc: &mut KappaConfig) -> Result<(), KappaError> {
    if let Err(e) = check_inputs(kpc) {
        match e {
            KappaError::SizeReorderRange => {
                eprintln!("Reordering size distribution range");
                swap(&mut kpc.amin, &mut kpc.amax);
            }
            KappaError::WavelengthReorderRange => {
                eprintln!("Reordering wavelength range");
                swap(&mut kpc.lmin, &mut kpc.lmax);
            }
            KappaError::SizeForceUnitSampling => {
                eprintln!("Setting na = 1 as amin = amax");
                kpc.na = 1;
            }
            KappaError::WavelengthForceUnitSampling => {
                eprintln!("Setting nlam = 1 as lmin = lmax");
                kpc.nlam = 1;
            }
            KappaError::SetSampling => {
                eprintln!("Setting sampling for size distribution");
                kpc.na = (((kpc.amax.log10() - kpc.amin.log10()) * 15.0 + 1.0) as usize).max(5);
            }
            KappaError::InvalidSizeParam(msg) => {
                return Err(KappaError::InvalidSizeParam(msg));
            }
            KappaError::ForceLogNormal => {
                eprintln!("Forcing log-normal distribution");
                kpc.sizedis = SizeDistribution::Lognorm;
            }
            KappaError::DisableSplit => {
                kpc.split = false;
            }
            e => return Err(e),
        }
    }
    Ok(())
}

fn check_inputs(kpc: &KappaConfig) -> Result<(), KappaError> {
    // Materials
    if kpc.nmat >= 20 {
        return Err(KappaError::TooManyMaterials);
    }
    if kpc.nmat == kpc.nmant && kpc.nmant > 0 {
        return Err(KappaError::NotEnoughMaterials);
    }

    // Porosity
    if kpc.pcore < 0.0 || kpc.pcore >= 1.0 || kpc.pmantle < 0.0 || kpc.pmantle >= 1.0 {
        return Err(KappaError::InvalidPorosity);
    }

    // Grain size distribution
    if kpc.amin <= 0.0 || kpc.amax <= 0.0 {
        return Err(KappaError::InvalidSizeInput);
    }
    if kpc.amin >= kpc.amax {
        // This error will actually be dealt with later
        return Err(KappaError::SizeReorderRange);
    }
    if kpc.na == 0 {
        // This error will actually be dealt with later
        return Err(KappaError::SamplingRequired);
    }
    if kpc.amin == kpc.amax && kpc.na != 1 {
        return Err(KappaError::SizeForceUnitSampling);
    }
    kpc.sizedis.validate(kpc)?;

    // Wavelength grid
    if kpc.lmin <= 0.0 || kpc.lmax <= 0.0 {
        return Err(KappaError::InvalidWavelengthInput);
    }
    if kpc.lmin > kpc.lmax {
        return Err(KappaError::WavelengthReorderRange);
    }
    if kpc.nlam <= 1 && kpc.lmin != kpc.lmax {
        return Err(KappaError::SamplingRequired);
    }
    if kpc.lmin == kpc.lmax && kpc.nlam != 1 {
        return Err(KappaError::WavelengthForceUnitSampling);
    }
    if kpc.nsparse > 0 {
        println!("Creating sparse scattering matrix file");
    }

    // DHS
    match kpc.method {
        KappaMethod::DHS => {
            if kpc.fmax < 0.0 || kpc.fmax >= 1.0 {
                return Err(KappaError::InvalidArgument);
            }
        }
        KappaMethod::MMF => {
            if kpc.mmf_struct > 3.0 {
                return Err(KappaError::InvalidArgument);
            }
            if kpc.mmf_struct <= 0.0 {
                return Err(KappaError::InvalidArgument);
            }
            if kpc.mmf_a0 >= kpc.amin {
                return Err(KappaError::InvalidArgument);
            }
        }
        KappaMethod::CDE => {
            if kpc.lmin <= 2.0 * PI * kpc.amax {
                eprintln!("WARNING: CDE requires Rayleigh limit!");
            }
        }
    }

    // Angular grid
    if kpc.nang % 2 == 1 {
        return Err(KappaError::InvalidArgument);
    }

    // Other
    if kpc.split && kpc.blend_only {
        eprintln!("WARNING: Turning off `-s` for `-blend_only`");
        return Err(KappaError::DisableSplit);
    }
    if kpc.split && kpc.sizedis != SizeDistribution::Apow {
        return Err(KappaError::InvalidArgument);
    }

    Ok(())
}

fn compute_kappa(kpc: &KappaConfig) -> Result<(), KappaError> {
    let (nf, ifmn) = if kpc.fmax == 0.0 { (1, 1) } else { (20, 12) };

    // bruggeman_blend();
    // maxwell_garnet_blend();
    todo!()
}

fn bruggeman_blend(abun: &[f64], e_in: &[Complex64]) -> Result<Complex64, KappaError> {
    let mut abunvac = 1.0 - abun.iter().sum::<f64>();
    let mvac = Complex::new(1.0, 0.0);
    if abunvac < 0.0 {
        if abunvac.abs() < 1e-5 {
            abunvac = 0.0;
        } else {
            return Err(KappaError::UnnormalizedAbundance);
        }
    }

    // Blend iteratively
    let mut mm = mvac;
    let mut tot: Complex64 = Complex::new(0.0, 0.0);
    let mut me: Complex64 = Complex::new(0.0, 0.0);
    for _ in 0..100 {
        tot = e_in
            .iter()
            .zip(abun.iter())
            .map(|(&e_in_j, &abun_j)| {
                (e_in_j.powi(2) - mm.powi(2)) / (e_in_j.powi(2) + 2.0 * mm.powi(2)) * abun_j
            })
            .sum::<Complex64>()
            + ((mvac.powi(2) - mm.powi(2)) / (mvac.powi(2) + 2.0 * mm.powi(2)) * abunvac);

        me = mm * ((2.0 * tot + 1.0) / (1.0 - tot)).sqrt();
        mm = me;
    }
    if tot.abs() / mm.abs() > 1e-6 {
        eprintln!("WARNING: Bruggeman blend not converged");
    }
    Ok(me)
}

fn maxwell_garnet_blend(m1: Complex64, m2: Complex64, vf_m: f64) -> (f64, f64) {
    let vf_c = 1.0 - vf_m;
    let me = m2.powi(2)
        * ((2.0 * m2.powi(2) + m1.powi(2) - 2.0 * vf_c * (m2.powi(2) - m1.powi(2)))
            / (2.0 * m2.powi(2) + m1.powi(2) + vf_c * (m2.powi(2) - m1.powi(2))));

    (me.sqrt().re, me.sqrt().im)
}
