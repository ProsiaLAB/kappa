//! The heart of `kappa`.
//!
//!

use std::f64::consts::PI;
use std::mem::swap;
use std::sync::Arc;

use anyhow::anyhow;
use anyhow::Result;
use ndarray::s;
use ndarray::Zip;
use num_complex::ComplexFloat;
use num_complex::{Complex, Complex64};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::io::write_sizedis_file;
use crate::types::RVecView;
use crate::types::{CVector, RMatrix, RVector, UVector};
use crate::utils::legendre::gauss_legendre;

#[derive(Debug)]
pub struct Material {
    pub key: String,
    pub kind: MaterialKind,
    pub n: f64,
    pub k: f64,
    pub rho: f64,
    pub mfrac: f64,
    pub re: RVector,
    pub im: RVector,
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
            re: RVector::zeros(nlam),
            im: RVector::zeros(nlam),
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
    pub ameans_file: [f64; 3],
    pub sizedis: SizeDistribution,
    pub wavelength_kind: WavelengthKind,
    pub lmin: f64,
    pub lmax: f64,
    pub nlam: usize,
    pub lam: RVector,
    pub iscatlam: UVector,
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
            wavelength_kind: WavelengthKind::CmdLine,
            lmin: 0.05,
            lmax: 1e4,
            nlam,
            lam: RVector::zeros(nlam),
            iscatlam: UVector::zeros(nlam),
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
    r: RVector,
    lam: RVecView<'a>,
    nlam: usize,
    mu: RVecView<'a>,
    e1_blend: RVecView<'a>,
    e2_blend: RVecView<'a>,
    p: Particle,
    nr: RVecView<'a>,
    method: KappaMethod,
    nf: usize,
    ifmn: usize,
    ns: usize,
    pcore: f64,
    rho_av: f64,
    wf: RVector,
    f: RVector,
    split: bool,
    nang: usize,
    chopangle: f64,
    qabsdqext_min: f64,
    xlim: f64,
    xlim_dhs: f64,
    mmf_a0: f64,
    mmf_struct: f64,
    mmf_kf: f64,
    mmfss: bool,
}

impl<'a> KappaState<'a> {
    pub fn new(
        kpc: &'a KappaConfig,
        ns: usize,
        nf: usize,
        r: RVector,
        nr: &'a RVector,
        mu: &'a RVector,
        e1_blend: &'a RVector,
        e2_blend: &'a RVector,
    ) -> Self {
        KappaState {
            r,
            lam: kpc.lam.view(),
            nlam: kpc.nlam,
            mu: mu.view(),
            e1_blend: e1_blend.view(),
            e2_blend: e2_blend.view(),
            p: Particle::default(),
            nr: nr.view(),
            method: KappaMethod::DHS,
            nf,
            ifmn: 0,
            ns,
            pcore: 0.0,
            rho_av: 0.0,
            wf: RVector::zeros(nf),
            f: RVector::zeros(nf),
            split: true,
            nang: kpc.nang,
            chopangle: 0.0,
            qabsdqext_min: 0.0,
            xlim: 1.0,
            xlim_dhs: 1.0,
            mmf_a0: 0.0,
            mmf_struct: 0.0,
            mmf_kf: 0.0,
            mmfss: false,
        }
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
        WavelengthKind::CmdLine => {
            kpc.lam = RVector::logspace(10.0, kpc.lmin, kpc.lmax, kpc.nlam);
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

    let mut r = RVector::zeros(ns);
    let mut nr = RVector::zeros(ns);
    let mut f = RVector::zeros(nf);
    let mut wf = RVector::zeros(nf);
    let mut e1mantle = RVector::zeros(kpc.nlam);
    let mut e2mantle = RVector::zeros(kpc.nlam);

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
    let mfrac_mantle = tot_mantle / tot;

    let aminlog = kpc.amin.log10();
    let amaxlog = kpc.amax.log10();
    let pow = -kpc.apow;

    if ns == 1 {
        // Just one size
        r[0] = 10.0_f64.powf((aminlog + amaxlog) / 2.0);
        nr[0] = r[0].powf(pow + 1.0); // should be 1/r[0]^3 ???  Not important.
    } else {
        let mut tot = 0.0;
        // Size distribution
        let nsf = (ns - 1) as f64;
        Zip::indexed(&mut r)
            .and(&mut nr)
            .for_each(|is, r_val, nr_val| {
                // Calculate size distribution values
                let isf = is as f64;
                *r_val = 10.0f64.powf(aminlog + (amaxlog - aminlog) * isf / nsf);

                // Apply size distribution and exponential calculation
                *nr_val = match kpc.sizedis {
                    // Log-normal or normal size distribution
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
                    // Power law size distribution
                    _ => r_val.powf(pow + 1.0),
                };

                // With the option `-d`m each computation is only a piece of the size grid
                if *r_val < kpc.amin || *r_val > kpc.amax {
                    *nr_val = 0.0;
                }

                // Volume normalization
                tot += *nr_val * r_val.powi(3);
            });

        // Normalize nr
        nr = 1.0 * nr / tot;
    }

    if kpc.write_grid {
        write_sizedis_file(kpc, ns, &r, &mut nr, tot)?;
    }

    // Copy the refractory index data for all materials into local arrays
    let mut e1 = RMatrix::zeros((kpc.nlam, kpc.nmat + 1));
    let mut e2 = RMatrix::zeros((kpc.nlam, kpc.nmat + 1));
    for im in 0..kpc.nmat {
        e1.slice_mut(s![.., im]).assign(&kpc.materials[im].re);
        e2.slice_mut(s![.., im]).assign(&kpc.materials[im].im);
    }

    let mut rho_mantle: f64 = 0.0;
    let mut vfrac_mantle = RVector::zeros(kpc.nmant);

    // Core: Turn mass fractions into volume fractions, compute rho_core
    let mtot_core = kpc
        .materials
        .iter()
        .filter(|m| m.kind == MaterialKind::Core)
        .fold(0.0, |sum, m| sum + m.mfrac);
    let mut vfrac_core = kpc
        .materials
        .iter()
        .filter(|m| m.kind == MaterialKind::Core)
        .map(|m| m.mfrac / m.rho)
        .collect::<RVector>();
    let vtot_core = vfrac_core.iter().sum::<f64>();
    let mut rho_core = mtot_core / vtot_core;
    vfrac_core /= vtot_core;
    if kpc.pcore > 0.0 {
        vfrac_core *= 1.0 - kpc.pcore;
        rho_core *= 1.0 - kpc.pcore;
    }

    // Mantle: Turn mass fractions to volume fractions, compute rho_mantle
    if kpc.nmant > 0 {
        let mtot_mantle = kpc
            .materials
            .iter()
            .filter(|m| m.kind == MaterialKind::Mantle)
            .fold(0.0, |sum, m| sum + m.mfrac);
        vfrac_mantle = kpc
            .materials
            .iter()
            .filter(|m| m.kind == MaterialKind::Mantle)
            .map(|m| m.mfrac / m.rho)
            .collect::<RVector>();
        let vtot_mantle = vfrac_mantle.iter().sum::<f64>();
        rho_mantle = mtot_mantle / vtot_mantle;
        vfrac_mantle /= vtot_mantle;
        if kpc.pmantle > 0.0 {
            vfrac_mantle *= 1.0 - kpc.pmantle;
            rho_mantle *= 1.0 - kpc.pmantle;
        }
    }

    let rho_av = if kpc.nmant == 0 {
        rho_core
    } else {
        rho_core / (1.0 + mfrac_mantle * (rho_core / rho_mantle - 1.0))
    };
    let vfrac_mantle_f = mfrac_mantle * rho_av / rho_mantle;

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
            let e_out = bruggeman_blend(&mfrac.to_vec(), &e_in.to_vec())?;
            e1_blend[il] = e_out.re;
            e2_blend[il] = e_out.im;
        }

        // Mantle
        if kpc.nmant > 0 {
            // We do have a mantle to add
            if (kpc.nmant == 1) && (kpc.pmantle == 0.0) {
                // No Blending needed inside the mantle - just copy e1 and e2
                // Since it is onyl one material, we know it is index nm
                e1_blend[il] = e1[[kpc.nmat, il]];
                e2_blend[il] = e2[[kpc.nmat, il]];
            } else {
                // Blend the mantle materials
                for im in 0..kpc.nmant {
                    e_in[im] = Complex::new(e1[[il, im + kpc.ncore]], e2[[il, im + kpc.ncore]]);
                }
                let e_out = bruggeman_blend(&vfrac_mantle.to_vec(), &e_in.to_vec())?;
                e1mantle[il] = e_out.re;
                e1mantle[il] = e_out.im;
            }
            let (e1mg, e2mg) = maxwell_garnet_blend(
                Complex::new(e1_blend[il], e2_blend[il]),
                Complex::new(e1mantle[il], e2mantle[il]),
                vfrac_mantle_f,
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
        return Ok(());
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
    for (j, val) in mu.iter_mut().enumerate() {
        let jf = j as f64;
        let theta = (jf - 0.5) / nangby2f * PI / 2.0;
        *val = theta.cos();
    }

    // Create shared memory for the variables
    let kps = Arc::new(KappaState::new(
        kpc, ns, nf, r, &nr, &mu, &e1_blend, &e2_blend,
    ));

    if kpc.split {
        let num_threads = rayon::current_num_threads();
        println!("Rayon thread pool size: {}", num_threads);
        (0..kpc.nlam).into_par_iter().for_each(|ilam| {
            let wvno = 2.0 * PI / kps.lam[ilam];
            let f11 = RVector::zeros(kpc.nang);
            let f12 = RVector::zeros(kpc.nang);
            let f22 = RVector::zeros(kpc.nang);
            let f33 = RVector::zeros(kpc.nang);
            let f34 = RVector::zeros(kpc.nang);
            let f44 = RVector::zeros(kpc.nang);
            let mut csca = 0.0;
            let mut cabs = 0.0;
            let mut cext = 0.0;
            let mut smat_nbad = 0;
            let mut mass = 0.0;
            let mut vol = 0.0;
            for is in 0..ns {
                let r1 = kps.r[is];
                let err = 0;
                let spheres = 0;
                let toolarge = 0;
                match kpc.method {
                    KappaMethod::DHS => {
                        let m_in = Complex::new(1.0, 0.0);
                        for ifmn in 0..nf {}
                    }
                    KappaMethod::MMF => {}
                    KappaMethod::CDE => {}
                }
            }
        });
    }

    // todo!()
    Ok(())
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
