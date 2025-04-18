//! The heart of `kappa`.
//! This module contains the main routines to compute opacities
//! and scattering matrices.

use std::cmp::Ordering;
use std::collections::HashSet;
use std::f64::consts::PI;
use std::mem::swap;
use std::process::exit;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering as AtomicOrdering;
use std::sync::Arc;

use anyhow::anyhow;
use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::s;
use num_complex::ComplexFloat;
use num_complex::{Complex, Complex64};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::components::get_lnk_data;
use crate::dhs::toon_ackerman_1981;
use crate::dhs::DHSConfig;
use crate::fractal::mean_scattering;
use crate::fractal::FractalConfig;
use crate::fractal::{FractalCutoff, FractalGeometry, FractalSolver};
use crate::io::{write_opacities, write_sizedis_file, write_wavelength_grid};
use crate::mie::de_rooij_1984;
use crate::mie::MieConfig;
use crate::types::{BVector, CVector, RMatrix, RVector};
use crate::utils::legendre::gauss_legendre;
use crate::utils::{prepare_sparse, regrid_lnk_data};

#[derive(Debug)]
pub struct Material {
    pub key: String,
    pub kind: MaterialKind,
    pub n: f64,
    pub k: f64,
    pub rho: f64,
    pub mfrac: f64,
    pub vfrac: f64,
    pub re: RVector,
    pub im: RVector,
    pub cmd: RefractiveIndexKind,
}

impl Default for Material {
    fn default() -> Self {
        let nlam = 300;
        Material {
            key: String::new(),
            kind: MaterialKind::Core,
            n: 0.0,
            k: 0.0,
            rho: 0.0,
            mfrac: 1.0,
            vfrac: 0.0,
            re: RVector::zeros(nlam),
            im: RVector::zeros(nlam),
            cmd: RefractiveIndexKind::Other,
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum MaterialKind {
    Core,
    Mantle,
}

#[derive(Debug, PartialEq)]
pub enum RefractiveIndexKind {
    File,
    CmdLine,
    Other,
}

#[derive(Debug)]
pub struct KappaConfig {
    pub amin: f64,
    pub amax: f64,
    pub apow: f64,
    pub amean: f64,
    pub asigma: f64,
    pub na: usize,
    pub ameans_file: [f64; 3],
    pub sizedis: SizeDistribution,
    pub wavelength_kind: WavelengthKind,
    pub lmin: f64,
    pub lmax: f64,
    pub nlam: usize,
    pub lam: RVector,
    pub sparse_indices: HashSet<usize>,
    pub nang: usize,
    pub chop_angle: f64,
    pub pcore: f64,
    pub pmantle: f64,
    pub nmat: usize,
    pub ncore: usize,
    pub nmant: usize,
    pub rho_core: f64,
    pub rho_mantle: f64,
    pub rho_av: f64,
    pub tot_mfrac_core: f64,
    pub tot_vfrac_core: f64,
    pub tot_mfrac_mantle: f64,
    pub tot_vfrac_mantle: f64,
    pub method: KappaMethod,
    pub fmax: f64,
    pub write_fits: bool,
    pub write_scatter: bool,
    pub write_grid: bool,
    pub for_radmc: bool,
    pub nsparse: usize,
    pub scatlammin: RVector,
    pub scatlammax: RVector,
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
    pub outdir: String,
}

impl Default for KappaConfig {
    fn default() -> Self {
        let nlam = 300;
        KappaConfig {
            amin: 0.05,
            amax: 3000.0,
            apow: 3.5,
            amean: 0.0,
            asigma: 0.0,
            na: 0,
            ameans_file: [0.0; 3],
            sizedis: SizeDistribution::Apow,
            wavelength_kind: WavelengthKind::Other,
            lmin: 0.05,
            lmax: 1e4,
            nlam,
            lam: RVector::zeros(nlam),
            sparse_indices: HashSet::new(),
            nang: 180,
            chop_angle: 0.0,
            pcore: 0.0,
            pmantle: 0.0,
            nmat: 0,
            ncore: 0,
            nmant: 0,
            rho_core: 0.0,
            rho_mantle: 0.0,
            rho_av: 0.0,
            tot_mfrac_core: 0.0,
            tot_vfrac_core: 0.0,
            tot_mfrac_mantle: 0.0,
            tot_vfrac_mantle: 0.0,
            method: KappaMethod::DHS,
            fmax: 0.8,
            write_fits: false,
            write_scatter: false,
            write_grid: false,
            for_radmc: false,
            nsparse: 0,
            scatlammin: RVector::zeros(30),
            scatlammax: RVector::zeros(30),
            nsubgrains: 5,
            mmf_struct: 0.0,
            mmf_a0: 0.0,
            mmf_kf: 0.0,
            mmf_ss: false,
            split: false,
            blend_only: false,
            xlim: 1e8,
            xlim_dhs: 1e4,
            materials: Vec::with_capacity(20),
            outdir: String::from("output"),
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
            key: "pyr-mg70".into(),
            kind: MaterialKind::Core,
            mfrac: 0.87,
            ..Default::default()
        });
        kpc.materials.push(Material {
            key: "c-z".into(),
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
            key: "astrosil".into(),
            kind: MaterialKind::Core,
            mfrac: 0.3291,
            ..Default::default()
        });
        kpc.materials.push(Material {
            key: "c-org".into(),
            kind: MaterialKind::Core,
            mfrac: 0.3966,
            ..Default::default()
        });
        kpc.materials.push(Material {
            key: "fes".into(),
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
            key: "h2o-w".into(),
            kind: MaterialKind::Core,
            mfrac: 0.2000,
            ..Default::default()
        });
        kpc
    }
}

struct KappaState<'a> {
    ns: usize,
    nf: usize,
    ifmn: usize,
    r: &'a RVector,
    nr: &'a RVector,
    f: &'a RVector,
    wf: &'a RVector,
    e1_blend: &'a RVector,
    e2_blend: &'a RVector,
    mu: &'a RVector,
}

#[derive(PartialEq, Debug)]
pub enum SizeDistribution {
    Apow,
    File,
    Normal,
    LogNormal,
}

#[derive(Debug)]
pub enum WavelengthKind {
    Other,
    File,
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

#[derive(Default)]
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
    pub trust: BVector,
    pub is_ok: bool,
    pub is_ok_lmin: f64,
}

struct KappaResult {
    mueller: Mueller,
    mass: f64,
    vol: f64,
    k_ext: f64,
    k_sca: f64,
    k_abs: f64,
    g: f64,
    trust: bool,
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
    /// Number of wavelengths.
    pub size: usize,
    /// Wavelengths.
    pub l0: RVector,
    /// Refractive index.
    pub n0: RVector,
    /// Extinction coefficient.
    pub k0: RVector,
}

/// Run the simulation.
pub fn run(kpc: &mut KappaConfig) -> Result<(), KappaError> {
    prepare_inputs(kpc)?;

    for material in &mut kpc.materials {
        match material.cmd {
            RefractiveIndexKind::CmdLine => {}
            RefractiveIndexKind::File | RefractiveIndexKind::Other => {
                let component = get_lnk_data(&material.key);
                (material.re, material.im) =
                    regrid_lnk_data(component.l0, component.n0, component.k0, &kpc.lam, true);
                material.rho = component.rho;
            }
        }
    }

    // Normalize the mass fractions
    (kpc.tot_mfrac_core, kpc.tot_mfrac_mantle) =
        kpc.materials
            .iter()
            .fold((0.0, 0.0), |(tot_core, tot_mantle), m| match m.kind {
                MaterialKind::Core => (tot_core + m.mfrac, tot_mantle),
                MaterialKind::Mantle => (tot_core, tot_mantle + m.mfrac),
            });
    kpc.materials.iter_mut().for_each(|m| match m.kind {
        MaterialKind::Core => {
            m.mfrac /= kpc.tot_mfrac_core;
        }
        MaterialKind::Mantle => {
            m.mfrac /= kpc.tot_mfrac_mantle;
        }
    });
    kpc.materials.iter_mut().for_each(|m| {
        m.vfrac = m.mfrac / m.rho;
    });
    (kpc.tot_vfrac_core, kpc.tot_vfrac_mantle) =
        kpc.materials
            .iter()
            .fold((0.0, 0.0), |(tot_core, tot_mantle), m| match m.kind {
                MaterialKind::Core => (tot_core + m.vfrac, tot_mantle),
                MaterialKind::Mantle => (tot_core, tot_mantle + m.vfrac),
            });

    // Turn mass fractions into volume fractions
    kpc.rho_core = kpc.tot_mfrac_core / kpc.tot_vfrac_core;
    kpc.rho_mantle = kpc.tot_mfrac_mantle / kpc.tot_vfrac_mantle;

    // Normalize volume fractions
    kpc.materials.iter_mut().for_each(|m| match m.kind {
        MaterialKind::Core => m.vfrac /= kpc.tot_vfrac_core,
        MaterialKind::Mantle => m.vfrac /= kpc.tot_vfrac_mantle,
    });

    if kpc.pcore > 0.0 {
        kpc.materials.iter_mut().for_each(|m| {
            if m.kind == MaterialKind::Core {
                m.vfrac *= 1.0 - kpc.pcore;
            }
        });
        kpc.rho_core *= 1.0 - kpc.pcore;
    }

    if kpc.pmantle > 0.0 {
        kpc.materials.iter_mut().for_each(|m| {
            if m.kind == MaterialKind::Mantle {
                m.vfrac *= 1.0 - kpc.pmantle;
            }
        });
        kpc.rho_mantle *= 1.0 - kpc.pmantle;
    }

    // Calculate average density of the whole grain
    if kpc.nmant == 0 {
        kpc.rho_av = kpc.rho_core;
    } else {
        kpc.rho_av =
            kpc.rho_core / (1.0 + kpc.tot_mfrac_mantle * (kpc.rho_core / kpc.rho_mantle - 1.0));
        kpc.tot_vfrac_mantle = kpc.tot_mfrac_mantle * kpc.rho_av / kpc.rho_mantle;
    }

    // Loop for splitting the output into files by grain size
    if kpc.split {
        let mut nsub = kpc.nsubgrains;
        if nsub % 2 == 0 {
            nsub += 1;
        }
        let afact = (kpc.amax / kpc.amin).powf(1.0 / kpc.na as f64);
        let afsub = afact.powf(1.0 / (nsub - 1) as f64);

        let bar = Arc::new(ProgressBar::new(kpc.nlam as u64));

        bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap()
                .progress_chars("█▓▒░ "),
        );

        (0..kpc.na).into_par_iter().try_for_each(|ia| {
            let iaf = ia as f64;
            let asplit = kpc.amin * afact.powf(iaf + 0.5);
            let nsubf = nsub as f64;
            let aminsplit = asplit * afsub.powf(-nsubf / 2.0);
            let amaxsplit = asplit * afsub.powf(nsubf / 2.0);
            let _ = compute_kappa(ia, aminsplit, amaxsplit, kpc)?;
            bar.inc(1); // Increment the progress bar inside the parallel loop
            Ok::<(), KappaError>(())
        })?;
        bar.finish_with_message("Processing complete!");
        exit(0);
    } else {
        let p = compute_kappa(0, kpc.amin, kpc.amax, kpc)?;
        if !p.is_ok {
            if kpc.mmf_ss {
                eprintln!("WARNING: opacities OK, but some F_nn,g_asym may not be accurate");
            } else {
                eprintln!("WARNING: opacities OK, but some F_nn,g_asym are set to zero");
            }
        }

        // TODO: FITS writer will be implemented later

        write_opacities(kpc, &p)?;
    }

    Ok(())
}

/// Pre-process the inputs collected from the command-line
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

    // Sort by `kind = MaterialKind::Core`
    kpc.materials.sort_by(|a, b| match (&a.kind, &b.kind) {
        (MaterialKind::Core, MaterialKind::Core) => Ordering::Equal,
        (MaterialKind::Core, _) => Ordering::Less,
        (_, MaterialKind::Core) => Ordering::Greater,
        _ => Ordering::Equal,
    });

    kpc.nmat = kpc.ncore + kpc.nmant;

    // Make logarithmic wavelength grid
    match kpc.wavelength_kind {
        WavelengthKind::Other => {
            kpc.lam = RVector::logspace(10.0, kpc.lmin.log10(), kpc.lmax.log10(), kpc.nlam);
        }
        WavelengthKind::File => {
            // already read from file
        }
    }

    // Prepare sparse scattering file
    if kpc.nsparse > 0 {
        prepare_sparse(kpc);
    }
    // Writing wavelength grid file
    if kpc.write_grid {
        write_wavelength_grid(kpc)?;
    }

    Ok(())
}

/// Calculate the opacities for a grain-size value.
fn compute_kappa(
    ia: usize,
    amin: f64,
    amax: f64,
    kpc: &KappaConfig,
) -> Result<Particle, KappaError> {
    let ns = kpc.na;
    let (nf, ifmn) = if kpc.fmax == 0.0 { (1, 1) } else { (20, 12) };

    let mut r = RVector::zeros(ns);
    let mut nr = RVector::zeros(ns);

    let mut e1mantle = RVector::zeros(kpc.nlam);
    let e2mantle = RVector::zeros(kpc.nlam);

    let aminlog = amin.log10();
    let amaxlog = amax.log10();
    let pow = -kpc.apow;

    let mut tot = 0.0;

    if ns == 1 {
        // Just one size
        r[0] = 10.0_f64.powf((aminlog + amaxlog) / 2.0);
        nr[0] = r[0].powf(pow + 1.0); // should be 1/r[0]^3 ???  Not important.
    } else {
        // Size distribution
        let nsf = (ns - 1) as f64;
        for (is, (r_val, nr_val)) in r.iter_mut().zip(nr.iter_mut()).enumerate() {
            let isf = is as f64;
            *r_val = 10.0f64.powf(aminlog + (amaxlog - aminlog) * isf / nsf);
            // Apply size distribution and exponential calculation
            *nr_val = match kpc.sizedis {
                // Power law size distribution
                SizeDistribution::Normal | SizeDistribution::LogNormal => {
                    if kpc.asigma < 0.0 {
                        let expo = 0.5 * ((*r_val - kpc.amean) / kpc.asigma).powf(2.0);
                        if expo > 99.0 {
                            0.0
                        } else {
                            (-expo).exp()
                        }
                    } else {
                        let expo = 0.5 * ((*r_val / kpc.amean).ln() / kpc.asigma).powf(2.0);
                        if expo > 99.0 {
                            0.0
                        } else {
                            (-expo).exp()
                        }
                    }
                }
                SizeDistribution::File => {
                    todo!() // TODO: Move this out of else block
                }
                _ => r_val.powf(pow + 1.0),
            };
            // With the option `-d`m each computation is only a piece of the size grid
            if *r_val < amin || *r_val > amax {
                *nr_val = 0.0;
            }
            // Volume normalization
            tot += *nr_val * r_val.powi(3);
        }
        // Normalize nr
        nr = 1.0 * nr / tot;
    }

    if kpc.write_grid && ia == 0 {
        write_sizedis_file(kpc, ns, &r, &mut nr, tot)?;
    }

    // Copy the refractory index data for all materials into local arrays
    let mut e1 = RMatrix::zeros((kpc.nlam, kpc.nmat));
    let mut e2 = RMatrix::zeros((kpc.nlam, kpc.nmat));
    for im in 0..kpc.nmat {
        e1.slice_mut(s![.., im]).assign(&kpc.materials[im].re);
        e2.slice_mut(s![.., im]).assign(&kpc.materials[im].im);
    }

    let mut e1_blend = RVector::zeros(kpc.nlam);
    let mut e2_blend = RVector::zeros(kpc.nlam);
    let mut e_in = CVector::zeros(kpc.ncore);

    // Mixing, for all wavelengths
    for il in 0..kpc.nlam {
        // Core
        if kpc.nmat == 1 && kpc.pcore == 0.0 {
            // Solid core, single material, nothing to blend for the core
            e1_blend[il] = e1[[0, il]];
            e2_blend[il] = e2[[0, il]];
        } else {
            // Blend the core materials
            for im in 0..kpc.ncore {
                e_in[im] = Complex::new(e1[[il, im]], e2[[il, im]]);
            }
            let e_out = bruggeman_blend(
                &kpc.materials
                    .iter()
                    .filter(|m| matches!(m.kind, MaterialKind::Core))
                    .map(|m| m.vfrac)
                    .collect::<Vec<_>>(),
                &e_in.to_vec(),
            )?;
            e1_blend[il] = e_out.re;
            e2_blend[il] = e_out.im;
        }

        // Mantle
        if kpc.nmant > 0 {
            // We do have a mantle to add
            if (kpc.nmant == 1) && (kpc.pmantle == 0.0) {
                // No Blending needed inside the mantle - just copy e1 and e2
                // Since it is only one material, we know it is index nm
                e1_blend[il] = e1[[kpc.nmat, il]];
                e2_blend[il] = e2[[kpc.nmat, il]];
            } else {
                // Blend the mantle materials
                for im in 0..kpc.nmant {
                    e_in[im] = Complex::new(e1[[il, im + kpc.ncore]], e2[[il, im + kpc.ncore]]);
                }
                let e_out = bruggeman_blend(
                    &kpc.materials
                        .iter()
                        .filter(|m| matches!(m.kind, MaterialKind::Mantle))
                        .map(|m| m.vfrac)
                        .collect::<Vec<_>>(),
                    &e_in.to_vec(),
                )?;
                e1mantle[il] = e_out.re;
                e1mantle[il] = e_out.im;
            }
            let (e1mg, e2mg) = maxwell_garnet_blend(
                Complex::new(e1_blend[il], e2_blend[il]), // Core
                Complex::new(e1mantle[il], e2mantle[il]), // Mantle
                kpc.tot_vfrac_mantle,
            );
            e1_blend[il] = e1mg;
            e2_blend[il] = e2mg;
        }
    }

    if kpc.blend_only {
        // Write .lnk file
        todo!()
    }
    // TODO: Now there is some `--print` command logic
    // TODO: that we will implement later

    if kpc.write_grid || kpc.blend_only {
        println!("Exiting after writing requested files");
        exit(0);
    }

    // Check how we are going to average over hollow sphere components
    let (f, wf) = if nf > 1 && kpc.fmax > 0.01 {
        // Get the weights for Gauss-Legendre integration
        gauss_legendre(0.01f64, kpc.fmax, nf)
    } else if kpc.fmax == 0.0 {
        (RVector::zeros(nf), RVector::ones(nf))
    } else {
        // Just a compact sphere, weight is 1
        (RVector::ones(nf) * kpc.fmax, RVector::ones(nf))
    };

    // Initialize mu
    let mut mu = RVector::zeros(kpc.nang);
    let nangby2f = (kpc.nang / 2) as f64;
    for (j, val) in mu.iter_mut().take(kpc.nang / 2).enumerate() {
        let jf = j as f64;
        let theta = (jf + 0.5) / nangby2f * PI / 2.0;
        *val = theta.cos();
    }

    // Create shared memory for the variables
    let kps = KappaState {
        ns,
        nf,
        ifmn,
        r: &r,
        nr: &nr,
        f: &f,
        wf: &wf,
        e1_blend: &e1_blend,
        e2_blend: &e2_blend,
        mu: &mu,
    };

    let p = {
        let bar = ProgressBar::new(kpc.nlam as u64);
        bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap()
                .progress_chars("█▓▒░ "),
        );
        let counter = Arc::new(AtomicUsize::new(0));
        let results: Result<Vec<KappaResult>, _> = if !kpc.split {
            let counter = Arc::clone(&counter);
            let bar = bar.clone();
            let res = (0..kpc.nlam)
                .into_par_iter()
                .map(|ilam| {
                    let r = over_wavelengths(ilam, &kps, kpc);
                    let prev = counter.fetch_add(1, AtomicOrdering::SeqCst);
                    bar.set_position((prev + 1) as u64);
                    r
                })
                .collect();
            bar.finish_with_message("Processing complete!");
            res
        } else {
            (0..kpc.nlam)
                .map(|ilam| over_wavelengths(ilam, &kps, kpc))
                .collect()
        };

        results.map(|lamrs| {
            let rho = lamrs[0].mass / lamrs[0].vol;
            let k_ext = lamrs.iter().map(|lamr| lamr.k_ext).collect::<RVector>();
            let k_sca = lamrs.iter().map(|lamr| lamr.k_sca).collect::<RVector>();
            let k_abs = lamrs.iter().map(|lamr| lamr.k_abs).collect::<RVector>();
            let g = lamrs.iter().map(|lamr| lamr.g).collect::<RVector>();
            let trust = lamrs.iter().map(|lamr| lamr.trust).collect::<BVector>();
            let mueller = lamrs
                .into_iter()
                .map(|lamr| lamr.mueller)
                .collect::<Vec<Mueller>>();

            let mut is_ok = true;
            let mut is_ok_lmin = 0.0;

            for ilam in 0..kpc.nlam {
                if !trust[ilam] {
                    is_ok = false;
                    is_ok_lmin = kpc.lam[ilam];
                };
            }

            Particle {
                rho,
                k_abs,
                k_ext,
                k_sca,
                g,
                f: mueller,
                trust,
                is_ok,
                is_ok_lmin,
                ..Default::default()
            }
        })?
    };

    Ok(p)
}

fn bruggeman_blend(abun: &[f64], e_in: &[Complex64]) -> Result<Complex64> {
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
        tot = ((mvac.powi(2) - mm.powi(2)) / (mvac.powi(2) + 2.0 * mm.powi(2))) * abunvac;
        tot += e_in
            .iter()
            .zip(abun.iter())
            .map(|(&e_in_j, &abun_j)| {
                ((e_in_j.powi(2) - mm.powi(2)) / (e_in_j.powi(2) + 2.0 * mm.powi(2))) * abun_j
            })
            .sum::<Complex64>();

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

/// Calculate the opacities for a wavelength value.
fn over_wavelengths(ilam: usize, kps: &KappaState, kpc: &KappaConfig) -> Result<KappaResult> {
    let qabsdqext_min = 1e-4;
    let wvno = 2.0 * PI / kpc.lam[ilam];
    let mut f11 = RVector::zeros(kpc.nang);
    let mut f12 = RVector::zeros(kpc.nang);
    let mut f22 = RVector::zeros(kpc.nang);
    let mut f33 = RVector::zeros(kpc.nang);
    let mut f34 = RVector::zeros(kpc.nang);
    let mut f44 = RVector::zeros(kpc.nang);
    let mut csca = 0.0;
    let mut cabs = 0.0;
    let mut cext = 0.0;
    let mut smat_nbad = 0;
    let mut mass = 0.0;
    let mut vol = 0.0;
    let nangf = kpc.nang as f64;
    let mut mie_f11 = RVector::zeros(kpc.nang);
    let mut mie_f12 = RVector::zeros(kpc.nang);
    let mut mie_f22 = RVector::zeros(kpc.nang);
    let mut mie_f33 = RVector::zeros(kpc.nang);
    let mut mie_f34 = RVector::zeros(kpc.nang);
    let mut mie_f44 = RVector::zeros(kpc.nang);
    for is in 0..kps.ns {
        let r1 = kps.r[is];
        let mut spheres = false;
        let mut too_large = false;
        match kpc.method {
            KappaMethod::DHS => {
                let m_in = Complex::new(1.0, 0.0);
                for ifn in 0..kps.nf {
                    let mut rad: f64;
                    if kps.f[ifn] == 0.0 {
                        spheres = true;
                    } else if r1 * wvno > kpc.xlim {
                        too_large = true;
                    } else {
                        rad = r1 / (1.0 - kps.f[ifn]).powf(1.0 / 3.0);
                        let rcore = rad * kps.f[ifn].powf(1.0 / 3.0);
                        let mconj = Complex::new(kps.e1_blend[ilam], -kps.e2_blend[ilam]);
                        let wvno_1 = wvno.min(kpc.xlim_dhs / rad);
                        let dhsc = DHSConfig {
                            r_core: rcore,
                            r_shell: rad,
                            wave_number: wvno_1,
                            r_indsh: mconj,
                            r_indco: m_in,
                            mu: kps.mu,
                            numang: kpc.nang / 2,
                            max_angle: kpc.nang,
                        };
                        let (csmie, mut cemie) = match toon_ackerman_1981(&dhsc) {
                            Ok(dhsr) => {
                                if spheres {
                                    rad = r1;
                                } else if too_large {
                                    rad = r1 / (1.0 - kps.f[kps.ifmn]).powf(1.0 / 3.0);
                                };
                                let cemie = dhsr.q_ext * PI * rad.powi(2);
                                let csmie = dhsr.q_sca * PI * rad.powi(2);
                                let factor = 2.0 * PI / csmie / wvno.powi(2);
                                for j in 0..(kpc.nang / 2) {
                                    mie_f11[j] = (dhsr.m1[[j, 0]] + dhsr.m0[[j, 0]]) * factor;
                                    mie_f12[j] = (dhsr.m1[[j, 0]] - dhsr.m0[[j, 0]]) * factor;
                                    mie_f22[j] = (dhsr.m1[[j, 0]] + dhsr.m0[[j, 0]]) * factor;
                                    mie_f33[j] = (dhsr.s10[[j, 0]]) * factor;
                                    mie_f34[j] = (-dhsr.d10[[j, 0]]) * factor;
                                    mie_f44[j] = (dhsr.s10[[j, 0]]) * factor;
                                    // Here we use the assumption that the grid is regular.  An adapted
                                    // grid is not possible if it is not symmetric around pi/2.
                                    mie_f11[kpc.nang - j - 1] =
                                        (dhsr.m1[[j, 1]] + dhsr.m0[[j, 1]]) * factor;
                                    mie_f12[kpc.nang - j - 1] =
                                        (dhsr.m1[[j, 1]] - dhsr.m0[[j, 1]]) * factor;
                                    mie_f22[kpc.nang - j - 1] =
                                        (dhsr.m1[[j, 1]] + dhsr.m0[[j, 1]]) * factor;
                                    mie_f33[kpc.nang - j - 1] = (dhsr.s10[[j, 1]]) * factor;
                                    mie_f34[kpc.nang - j - 1] = (-dhsr.d10[[j, 1]]) * factor;
                                    mie_f44[kpc.nang - j - 1] = (dhsr.s10[[j, 1]]) * factor;
                                }
                                // exit(0);
                                (csmie, cemie)
                            }
                            Err(e) => {
                                println!("DHS error: {:?}", e);
                                rad = r1;
                                let rmie = rad;
                                let lmie = kpc.lam[ilam];
                                let e1_mie = kps.e1_blend[ilam];
                                let e2_mie = kps.e2_blend[ilam];

                                let thmin = 180.0 * (1.0 - 0.5) / nangf;
                                let thmax = 180.0 * (nangf - 0.5) / nangf;
                                let radius = if rmie / lmie < 5000.0 {
                                    rmie
                                } else {
                                    rmie / 5000.0
                                };
                                let miec = MieConfig {
                                    nangle: kpc.nang,
                                    delta: 1e-8,
                                    thmin,
                                    thmax,
                                    step: (thmax - thmin) / (nangf - 1.0),
                                    lam: lmie,
                                    cmm: Complex::new(e1_mie, e2_mie),
                                    rad: radius,
                                };
                                let mier = de_rooij_1984(&miec)?;
                                let cemie = mier.c_ext;
                                let csmie = mier.c_sca;
                                (csmie, cemie)
                            }
                        };

                        let mut tot = 0.0;
                        let mut tot2 = 0.0;
                        for j in 0..kpc.nang {
                            let jf = j as f64;
                            tot += mie_f11[j] * (PI * (jf + 0.5) / nangf).sin();
                            tot2 += (PI * (jf + 0.5) / nangf).sin();
                        }
                        mie_f11[0] += (tot2 - tot) / (PI * 0.5 / nangf).sin();
                        if mie_f11[0] < 0.0 {
                            mie_f11[0] = 0.0;
                        }
                        for j in 0..kpc.nang {
                            f11[j] += kps.wf[ifn] * kps.nr[is] * mie_f11[j] * csmie;
                            f12[j] += kps.wf[ifn] * kps.nr[is] * mie_f12[j] * csmie;
                            f22[j] += kps.wf[ifn] * kps.nr[is] * mie_f22[j] * csmie;
                            f33[j] += kps.wf[ifn] * kps.nr[is] * mie_f33[j] * csmie;
                            f34[j] += kps.wf[ifn] * kps.nr[is] * mie_f34[j] * csmie;
                            f44[j] += kps.wf[ifn] * kps.nr[is] * mie_f44[j] * csmie;
                        }
                        let mut camie = cemie - csmie;
                        if camie < cemie * qabsdqext_min {
                            camie *= 1e-4;
                            cemie = camie + csmie;
                        }
                        cext += kps.wf[ifn] * kps.nr[is] * cemie;
                        csca += kps.wf[ifn] * kps.nr[is] * csmie;
                        cabs += kps.wf[ifn] * kps.nr[is] * camie;
                        mass += kps.wf[ifn] * kps.nr[is] * kpc.rho_av * 4.0 * PI * r1.powi(3) / 3.0;
                        vol += kps.wf[ifn] * kps.nr[is] * 4.0 * PI * r1.powi(3) / 3.0;
                    }
                }
            }
            KappaMethod::MMF => {
                let m_mono = 4.0 * PI / 3.0 * kpc.mmf_a0.powi(3) * kpc.rho_av;
                let v_agg = 4.0 * PI / 3.0 * r1.powi(3);
                let m_agg = v_agg * kpc.rho_av;
                let nmono = m_agg / m_mono;
                let dfrac = if kpc.mmf_struct > 1.0 {
                    kpc.mmf_struct
                } else {
                    3.0 * nmono.ln() / (nmono - kpc.mmf_struct).ln()
                };
                let kfrac = if kpc.mmf_kf > 0.0 {
                    kpc.mmf_kf
                } else {
                    (5.0 / 3.0).powf(dfrac)
                };
                let fracc = FractalConfig {
                    solver: FractalSolver::ModifiedMeanField,
                    cutoff: FractalCutoff::Gaussian,
                    geometry: FractalGeometry::Tazaki,
                    nang: kpc.nang / 2 + 1,
                    pn: nmono,
                    r0: kpc.mmf_a0,
                    df: dfrac,
                    k0: kfrac,
                    lmd: kpc.lam[ilam],
                    refrel: Complex::new(kps.e1_blend[ilam], kps.e2_blend[ilam]),
                    ..Default::default()
                };
                let fracr = mean_scattering(&fracc)?;
                if fracr.dphi > 1.0 {
                    smat_nbad += 1;
                }
                let factor = 4.0 * PI / wvno.powi(2);
                for j in 0..kpc.nang {
                    f11[j] +=
                        kps.nr[is] * 0.5 * (fracr.smat[[0, j]] + fracr.smat[[0, j + 1]]) * factor;
                    f12[j] +=
                        kps.nr[is] * 0.5 * (fracr.smat[[1, j]] + fracr.smat[[1, j + 1]]) * factor;
                    f22[j] +=
                        kps.nr[is] * 0.5 * (fracr.smat[[0, j]] + fracr.smat[[0, j + 1]]) * factor;
                    f33[j] +=
                        kps.nr[is] * 0.5 * (fracr.smat[[2, j]] + fracr.smat[[2, j + 1]]) * factor;
                    f34[j] +=
                        kps.nr[is] * 0.5 * (fracr.smat[[3, j]] + fracr.smat[[3, j + 1]]) * factor;
                    f44[j] +=
                        kps.nr[is] * 0.5 * (fracr.smat[[2, j]] + fracr.smat[[2, j + 1]]) * factor;
                }
                cext += kps.nr[is] * fracr.c_ext;
                csca += kps.nr[is] * fracr.c_sca;
                cabs += kps.nr[is] * fracr.c_abs;
                mass += kps.nr[is] * m_agg;
                vol += kps.nr[is] * v_agg;
            }
            KappaMethod::CDE => {
                let v = 4.0 * PI * r1.powi(3) / 3.0;
                vol += kps.nr[is] * v;
                mass += kps.nr[is] * kpc.rho_av * v;
                let m = Complex::new(kps.e1_blend[ilam], kps.e2_blend[ilam]);
                let dcabs = 2.0 * wvno * v * (m.powi(2) / (m.powi(2) - 1.0) * m.powi(2).ln()).im;
                let dcsca = if (kps.e2_blend[ilam] / kps.e1_blend[ilam]).abs() < 1e-6 {
                    // Non-absorbing case
                    let mre = m.re;
                    wvno.powi(4) * v.powi(2) / (3.0 * PI) * (mre.powi(2) - 1.0 - mre.powi(2).ln())
                } else {
                    // Absorbing case
                    wvno.powi(4) * v.powi(2) * (m.powi(2) - 1.0).abs().powi(2)
                        / (3.0 * PI * m.powi(2).im)
                        * (m.powi(2) / (m.powi(2) - 1.0) * m.powi(2).ln()).im
                };
                cabs += kps.nr[is] * dcabs;
                csca += kps.nr[is] * dcsca;
                cext += kps.nr[is] * (dcabs + dcsca);
                if is == kps.ns - 1 {
                    // Final size, cscat is complete at this point
                    // Compute the scattering matrix deep in the Rayleigh limit
                    let lmie = kpc.lam[ilam];
                    let rmie = lmie / 1e3;
                    let e1_mie = kps.e1_blend[kpc.nlam - 1];
                    let e2_mie = kps.e2_blend[kpc.nlam - 1];
                    let thmin = 180.0 * (1.0 - 0.5) / nangf;
                    let thmax = 180.0 * (nangf - 0.5) / nangf;
                    let miec = MieConfig {
                        nangle: kpc.nang,
                        delta: 1e-8,
                        thmin,
                        thmax,
                        step: (thmax - thmin) / (nangf - 1.0),
                        lam: lmie,
                        cmm: Complex::new(e1_mie, e2_mie),
                        rad: rmie,
                    };
                    let mier = de_rooij_1984(&miec)?;
                    let mie_f22 = mier.f_11.clone();
                    let mie_f44 = mier.f_33.clone();
                    for j in 0..kpc.nang {
                        f11[j] = mier.f_11[j] * csca;
                        f12[j] = mier.f_12[j] * csca;
                        f22[j] = mie_f22[j] * csca;
                        f33[j] = mier.f_33[j] * csca;
                        f34[j] = mier.f_34[j] * csca;
                        f44[j] = mie_f44[j] * csca;
                    }
                }
            }
        }
    }

    // if ilam == 0 {
    //     let rho = mass / vol;
    // } // TODO: This can be moved out of parallel loop
    let mut k_ext = 1e4 * cext / mass;
    let mut k_sca = 1e4 * csca / mass;
    let k_abs = 1e4 * cabs / mass;
    f11 /= csca;
    f12 /= csca;
    f22 /= csca;
    f33 /= csca;
    f34 /= csca;
    f44 /= csca;

    if kpc.chop_angle > 0.0 {
        // Flatten the peak
        let ichop = (kpc.chop_angle / (180.0 / nangf)) as usize;
        if ichop > 1 {
            for i in 0..ichop {
                f11[i] = f11[ichop + 1];
                f12[i] = f12[ichop + 1];
                f22[i] = f22[ichop + 1];
                f33[i] = f33[ichop + 1];
                f34[i] = f34[ichop + 1];
                f44[i] = f44[ichop + 1];
            }
        }
        let mut tot = 0.0;
        let mut tot2 = 0.0;
        for j in 0..kpc.nang {
            let jf = j as f64;
            tot += f11[j] * (PI * (jf + 0.5) / nangf).sin() * (PI / nangf) * 2.0 * PI;
            tot2 += (PI * (jf + 0.5) / nangf).sin() * (PI / nangf) * 2.0 * PI;
        }
        for j in 0..kpc.nang {
            f11[j] = f11[j] * tot2 / tot;
            f12[j] = f12[j] * tot2 / tot;
            f22[j] = f22[j] * tot2 / tot;
            f33[j] = f33[j] * tot2 / tot;
            f34[j] = f34[j] * tot2 / tot;
            f44[j] = f44[j] * tot2 / tot;
        }
        k_sca *= tot / tot2;
        k_ext = k_sca + k_abs;
    }
    let mut tot = 0.0;
    let mut g = 0.0;
    let mut test_scat = true;

    if smat_nbad > 0 {
        test_scat = false;
    }
    for i in 0..kpc.nang {
        let ir = i as f64;
        g += f11[i] * (PI * (ir + 0.5) / nangf).cos() * (PI * (ir + 0.5) / nangf).sin();
        tot += f11[i] * (PI * (ir + 0.5) / nangf).sin();
    }
    g /= tot;
    if smat_nbad > 0 && !kpc.mmf_ss {
        g = 0.0;
        for i in 0..kpc.nang {
            f11[i] = 0.0;
            f12[i] = 0.0;
            f22[i] = 0.0;
            f33[i] = 0.0;
            f34[i] = 0.0;
            f44[i] = 0.0;
        }
    }

    let mueller = Mueller {
        f11,
        f12,
        f22,
        f33,
        f34,
        f44,
    };

    let lamr = KappaResult {
        mueller,
        mass,
        vol,
        k_ext,
        k_sca,
        k_abs,
        g,
        trust: test_scat,
    };

    Ok(lamr)
}
