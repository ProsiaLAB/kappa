//! The heart of `kappa`.
//!
//!

use std::f64::consts::PI;
use std::mem::swap;

use anyhow::anyhow;
use anyhow::Result;
use ndarray::Array1;
use num_complex::ComplexFloat;
use num_complex::{Complex, Complex64};

use crate::types::RVector;

#[derive(Debug)]
pub struct Material {
    pub key: String,
    pub kind: MaterialKind,
    pub cmd: bool,
    pub n: f64,
    pub k: f64,
    pub rho: f64,
    pub mfrac: f64,
}

impl Default for Material {
    fn default() -> Self {
        Material {
            key: String::new(),
            kind: MaterialKind::Core,
            cmd: false,
            n: 0.0,
            k: 0.0,
            rho: 0.0,
            mfrac: 1.0,
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum MaterialKind {
    Core,
    Mantle,
}

#[derive(Debug)]
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
    pub ncore: usize,
    pub nmant: usize,
    pub method: KappaMethod,
    pub fmax: f64,
    pub write_fits: bool,
    pub write_scatter: bool,
    pub write_grid: bool,
    pub for_radmc: bool,
    pub nsparse: usize,
    pub nsubgrains: usize,
    pub mmf_struct: f64,
    pub mmf_a0: f64,
    pub mmf_kf: f64,
    pub mmf_ss: bool,
    pub split: bool,
    pub blend_only: bool,
    pub xlim: f64,
    pub xlim_dhs: f64,
    pub materials: Vec<Material>,
    pub outdir: Option<String>,
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
            ncore: 0,
            nmant: 0,
            method: KappaMethod::DHS,
            fmax: 0.8,
            write_fits: false,
            write_scatter: false,
            write_grid: false,
            for_radmc: false,
            nsparse: 0,
            nsubgrains: 0,
            mmf_struct: 0.0,
            mmf_a0: 0.0,
            mmf_kf: 0.0,
            mmf_ss: false,
            split: true,
            blend_only: false,
            xlim: 1.0,
            xlim_dhs: 1.0,
            materials: Vec::with_capacity(20),
            outdir: Option::Some(String::from("output")),
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
        let mut kpc = KappaConfig {
            nmat: 2,
            ncore: 2,
            pcore: 0.25,
            ..Default::default()
        };
        kpc.materials.push(Material {
            key: "pyr-mg70".to_string(),
            kind: MaterialKind::Core,
            mfrac: 0.87,
            ..Default::default()
        });
        kpc.materials.push(Material {
            key: "c-z".to_string(),
            kind: MaterialKind::Core,
            mfrac: 0.13,
            ..Default::default()
        });
        kpc
    }

    fn dsharp_no_ice() -> Self {
        let mut kpc = KappaConfig {
            nmat: 3,
            ncore: 3,
            ..Default::default()
        };
        kpc.materials.push(Material {
            key: "astrosil".to_string(),
            kind: MaterialKind::Core,
            mfrac: 0.3291,
            ..Default::default()
        });
        kpc.materials.push(Material {
            key: "c-org".to_string(),
            kind: MaterialKind::Core,
            mfrac: 0.3966,
            ..Default::default()
        });
        kpc.materials.push(Material {
            key: "fes".to_string(),
            kind: MaterialKind::Core,
            mfrac: 0.0743,
            ..Default::default()
        });
        kpc
    }

    fn dsharp() -> Self {
        let mut kpc = KappaConfig::dsharp_no_ice();
        kpc.nmat += 1;
        kpc.ncore += 1;
        kpc.materials.push(Material {
            key: "h2o-w".to_string(),
            kind: MaterialKind::Core,
            mfrac: 0.2000,
            ..Default::default()
        });
        kpc
    }
}

#[derive(PartialEq, Debug)]
pub enum SizeDistribution {
    Apow,
    File,
    Normal,
    LogNormal,
}

impl SizeDistribution {
    fn validate(&self, kpc: &KappaConfig) -> Result<(), KappaError> {
        match self {
            SizeDistribution::Apow if kpc.apow < 0.0 => Err(anyhow!("UnexpectedApow").into()),
            SizeDistribution::Normal => {
                if kpc.amean <= 0.0 {
                    return Err(
                        anyhow!("`amean` must be positive for (log-)normal distribution").into(),
                    );
                }
                if kpc.asigma == 0.0 {
                    return Err(
                        anyhow!("`asigma` cannot be zero for (log-)normal distribution").into(),
                    );
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

#[derive(Debug)]
pub enum KappaMethod {
    DHS,
    MMF,
    CDE,
}

#[derive(Debug)]
pub enum KappaError {
    ZeroMassFraction,
    ForceLogNormal,
    Other(String),
}

impl From<anyhow::Error> for KappaError {
    fn from(err: anyhow::Error) -> Self {
        KappaError::Other(err.to_string())
    }
}

/// Mueller matrix structure
pub struct Mueller {
    pub f11: RVector,
    pub f12: RVector,
    pub f22: RVector,
    pub f33: RVector,
    pub f44: RVector,
    pub f34: RVector,
}

/// Particle
pub struct Particle {
    pub rv: f64,
    pub rvmin: f64,
    pub rvmax: f64,
    pub rho: f64,
    pub k_abs: RVector,
    pub k_sca: RVector,
    pub k_ext: RVector,
    pub g: RVector,
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
    pub l0: RVector,
    /// Refractive index.
    pub n0: RVector,
    /// Extinction coefficient.
    pub k0: RVector,
}

pub fn run(kpc: &mut KappaConfig) -> Result<()> {
    prepare_inputs(kpc)?;
    // println!("{:?}", kpc);
    compute_kappa(kpc)?;
    Ok(())
}

fn prepare_inputs(kpc: &mut KappaConfig) -> Result<()> {
    // Materials
    if kpc.nmat >= 20 {
        return Err(anyhow!("TooManyMaterials"));
    }
    if kpc.nmat == kpc.nmant && kpc.nmant > 0 {
        return Err(anyhow!("NoCoreMaterial"));
    }

    // Porosity
    if kpc.pcore < 0.0 || kpc.pcore >= 1.0 || kpc.pmantle < 0.0 || kpc.pmantle >= 1.0 {
        return Err(anyhow!("InvalidPorosity"));
    }

    // Grain size distribution
    if kpc.amin <= 0.0 || kpc.amax <= 0.0 {
        return Err(anyhow!("InvalidSizeInput"));
    }
    if kpc.amin >= kpc.amax {
        swap(&mut kpc.amin, &mut kpc.amax);
    }
    if kpc.na == 0 {
        kpc.na = (((kpc.amax.log10() - kpc.amin.log10()) * 15.0 + 1.0) as usize).max(5);
    }
    if kpc.amin == kpc.amax && kpc.na != 1 {
        kpc.na = 1;
    }
    match kpc.sizedis.validate(kpc) {
        Ok(()) => {}
        Err(KappaError::ForceLogNormal) => {
            kpc.sizedis = SizeDistribution::LogNormal;
        }
        Err(_) => return Err(anyhow!("InvalidSizeParam")),
    }

    // Wavelength grid
    if kpc.lmin <= 0.0 || kpc.lmax <= 0.0 {
        return Err(anyhow!("InvalidWavelengthInput"));
    }
    if kpc.lmin > kpc.lmax {
        swap(&mut kpc.lmin, &mut kpc.lmax);
    }
    if kpc.nlam <= 1 && kpc.lmin != kpc.lmax {
        return Err(anyhow!("SamplingRequired"));
    }
    if kpc.lmin == kpc.lmax && kpc.nlam != 1 {
        kpc.nlam = 1
    }
    if kpc.nsparse > 0 {
        println!("Creating sparse scattering matrix file");
    }

    // DHS
    match kpc.method {
        KappaMethod::DHS => {
            if kpc.fmax < 0.0 || kpc.fmax >= 1.0 {
                return Err(anyhow!("InvalidArgument: fmax"));
            }
        }
        KappaMethod::MMF => {
            if kpc.mmf_struct > 3.0 {
                return Err(anyhow!("`mmf_struct` must be between 1 and 3"));
            }
            if kpc.mmf_struct <= 0.0 {
                return Err(anyhow!("`mmf_struct` must be positive"));
            }
            if kpc.mmf_a0 >= kpc.amin {
                return Err(anyhow!("`mmf_a0` must be smaller than `amin`"));
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
        return Err(anyhow!("`nang` must be even"));
    }

    // Other
    if kpc.split && kpc.blend_only {
        eprintln!("WARNING: Turning off `-s` for `-blend_only`");
        kpc.split = false;
    }
    if kpc.split && kpc.sizedis != SizeDistribution::Apow {
        return Err(anyhow!(
            "`split` is only supported for powerlaw size distribution"
        ));
    }

    Ok(())
}

fn compute_kappa(kpc: &KappaConfig) -> Result<()> {
    let ns = kpc.na;
    let (nf, ifmn) = if kpc.fmax == 0.0 { (1, 1) } else { (20, 12) };

    let mut r: RVector = Array1::zeros(ns);
    let mut nr: RVector = Array1::zeros(ns);
    let mut f: RVector = Array1::zeros(nf);
    let mut wf: RVector = Array1::zeros(nf);
    let mut e1mantle: RVector = Array1::zeros(kpc.nlam);
    let mut e2mantle: RVector = Array1::zeros(kpc.nlam);

    // Normalize the mass fractions
    let (tot, tot_mantle) = kpc
        .materials
        .iter()
        .fold((0.0, 0.0), |(sum, sum_mantle), m| {
            let sum_mantle = if m.kind == MaterialKind::Mantle {
                sum_mantle + m.mfrac
            } else {
                sum_mantle
            };
            (sum + m.mfrac, sum_mantle)
        });

    let mfrac: RVector = kpc.materials.iter().map(|m| m.mfrac / tot).collect();
    let mfrac_mantle: RVector = kpc
        .materials
        .iter()
        .filter(|m| m.kind == MaterialKind::Mantle)
        .map(|m| m.mfrac / tot_mantle)
        .collect();

    let aminlog = kpc.amin.log10();
    let amaxlog = kpc.amax.log10();
    let pow = -kpc.apow;

    if ns == 1 {
        // Just one size
        r[0] = 10.0_f64.powf((aminlog + amaxlog) / 2.0);
        nr[0] = r[0].powf(pow + 1.0); // should be 1/r[0]^3 ???  Not important.
    }

    // bruggeman_blend();
    // maxwell_garnet_blend();
    todo!()
}

pub fn bruggeman_blend(abun: &[f64], e_in: &[Complex64]) -> Result<Complex64> {
    let mut abunvac = 1.0 - abun.iter().sum::<f64>();
    let mvac = Complex::new(1.0, 0.0);
    if abunvac < 0.0 {
        if abunvac.abs() < 1e-5 {
            abunvac = 0.0;
        } else {
            return Err(anyhow!("UnnormalizedAbundance"));
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

pub fn maxwell_garnet_blend(m1: Complex64, m2: Complex64, vf_m: f64) -> (f64, f64) {
    let vf_c = 1.0 - vf_m;
    let me = m2.powi(2)
        * ((2.0 * m2.powi(2) + m1.powi(2) - 2.0 * vf_c * (m2.powi(2) - m1.powi(2)))
            / (2.0 * m2.powi(2) + m1.powi(2) + vf_c * (m2.powi(2) - m1.powi(2))));

    (me.sqrt().re, me.sqrt().im)
}
