//! The heart of `kappa`.
//! This module contains the main routines to compute opacities
//! and scattering matrices.

use std::cmp::Ordering;
use std::collections::HashSet;
use std::f64::consts::PI;
use std::fs;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::mem::swap;
use std::path::Path;
use std::process::exit;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering as AtomicOrdering;

use anyhow::Result;
use anyhow::anyhow;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::s;
use num_complex::ComplexFloat;
use num_complex::{Complex, Complex64};
use prosia_extensions::types::{BVector, CVector, RMatrix, RVector};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::components::get_lnk_data;
use crate::dhs::DHSConfig;
use crate::dhs::toon_ackerman_1981;
use crate::fractal::FractalConfig;
use crate::fractal::mean_scattering;
use crate::fractal::{FractalCutoff, FractalGeometry, FractalSolver};
use crate::io::read_lnk_file;
use crate::mie::MieConfig;
use crate::mie::de_rooij_1984;
use crate::utils::legendre::gauss_legendre;
use crate::utils::regrid_lnk_data;

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

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize, Clone))]
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

#[allow(clippy::struct_excessive_bools)]
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
    pub rho_core: f64,
    pub rho_mantle: f64,
    pub rho_av: f64,
    pub total_mfrac: f64,
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
            sizedis: SizeDistribution::PowerLaw,
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
            rho_core: 0.0,
            rho_mantle: 0.0,
            rho_av: 0.0,
            total_mfrac: 0.0,
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

impl KappaConfig {
    #[must_use]
    pub fn ncore(&self) -> usize {
        self.materials
            .iter()
            .filter(|m| matches!(m.kind, MaterialKind::Core))
            .count()
    }

    #[must_use]
    pub fn nmant(&self) -> usize {
        self.materials
            .iter()
            .filter(|m| matches!(m.kind, MaterialKind::Mantle))
            .count()
    }

    /// Pre-process the inputs collected from the command-line
    ///
    /// # Errors
    /// Returns a [`KappaError`] if there is an issue with parsing the arguments,
    /// such as missing required arguments, invalid argument formats, or unknown options.
    pub fn prepare_inputs(&mut self) -> Result<()> {
        // Materials
        if self.materials.len() >= 20 {
            return Err(anyhow!("TooManyMaterials"));
        }
        if self.materials.len() == self.nmant() && self.nmant() > 0 {
            return Err(anyhow!("NoCoreMaterial"));
        }

        // Porosity
        if self.pcore < 0.0 || self.pcore >= 1.0 || self.pmantle < 0.0 || self.pmantle >= 1.0 {
            return Err(anyhow!("InvalidPorosity"));
        }

        // Grain size distribution
        if self.amin <= 0.0 || self.amax <= 0.0 {
            return Err(anyhow!("InvalidSizeInput"));
        }
        if self.amin >= self.amax {
            swap(&mut self.amin, &mut self.amax);
        }
        if self.na == 0 {
            self.na = (((self.amax.log10() - self.amin.log10()) * 15.0 + 1.0) as usize).max(5);
        }
        if (self.amin - self.amax).abs() < f64::EPSILON && self.na != 1 {
            self.na = 1;
        }
        match self.validate_size_distribution() {
            Ok(()) => {}
            Err(KappaError::ForceLogNormal) => {
                self.sizedis = SizeDistribution::LogNormal;
            }
            Err(_) => return Err(anyhow!("InvalidSizeParam")),
        }

        // Wavelength grid
        if self.lmin <= 0.0 || self.lmax <= 0.0 {
            return Err(anyhow!("InvalidWavelengthInput"));
        }
        if self.lmin > self.lmax {
            swap(&mut self.lmin, &mut self.lmax);
        }
        if self.nlam <= 1 && (self.lmin - self.lmax).abs() > f64::EPSILON {
            return Err(anyhow!("SamplingRequired"));
        }
        if (self.lmin - self.lmax).abs() < f64::EPSILON && self.nlam != 1 {
            self.nlam = 1;
        }
        if self.nsparse > 0 {
            println!("Creating sparse scattering matrix file");
        }

        // DHS
        match self.method {
            KappaMethod::DHS => {
                if self.fmax < 0.0 || self.fmax >= 1.0 {
                    return Err(anyhow!("InvalidArgument: fmax"));
                }
            }
            KappaMethod::MMF => {
                if self.mmf_struct > 3.0 {
                    return Err(anyhow!("`mmf_struct` must be between 1 and 3"));
                }
                if self.mmf_struct <= 0.0 {
                    return Err(anyhow!("`mmf_struct` must be positive"));
                }
                if self.mmf_a0 >= self.amin {
                    return Err(anyhow!("`mmf_a0` must be smaller than `amin`"));
                }
            }
            KappaMethod::CDE => {
                if self.lmin <= 2.0 * PI * self.amax {
                    eprintln!("WARNING: CDE requires Rayleigh limit!");
                }
            }
        }

        // Angular grid
        if self.nang % 2 == 1 {
            return Err(anyhow!("`nang` must be even"));
        }

        // Other
        if self.split && self.blend_only {
            eprintln!("WARNING: Turning off `-s` for `-blend_only`");
            self.split = false;
        }
        if self.split && self.sizedis != SizeDistribution::PowerLaw {
            return Err(anyhow!(
                "`split` is only supported for powerlaw size distribution"
            ));
        }

        // Sort by `kind = MaterialKind::Core`
        self.materials.sort_by(|a, b| match (&a.kind, &b.kind) {
            (MaterialKind::Core, MaterialKind::Core) => Ordering::Equal,
            (MaterialKind::Core, _) => Ordering::Less,
            (_, MaterialKind::Core) => Ordering::Greater,
            _ => Ordering::Equal,
        });

        // Make logarithmic wavelength grid
        match self.wavelength_kind {
            WavelengthKind::Other => {
                self.lam = RVector::logspace(10.0, self.lmin.log10(), self.lmax.log10(), self.nlam);
            }
            WavelengthKind::File => {
                // already read from file
            }
        }

        // Prepare sparse scattering file
        if self.nsparse > 0 {
            self.prepare_sparse();
        }
        // Writing wavelength grid file
        if self.write_grid {
            self.write_wavelength_grid()?;
        }

        Ok(())
    }

    /// Initialize the `KappaConfig` struct.
    ///
    /// # Errors
    /// Returns a [`KappaError`] if there is an issue with parsing the arguments,
    /// such as missing required arguments, invalid argument formats, or unknown options.
    pub fn initialize(&mut self) -> Result<(), KappaError> {
        for material in &mut self.materials {
            match material.cmd {
                RefractiveIndexKind::CmdLine => {}
                RefractiveIndexKind::File => {
                    let component = read_lnk_file(&material.key, None)?;
                    let l0_slice: &[f64] = &component.l0.to_vec();
                    let n0_slice: &[f64] = &component.n0.to_vec();
                    let k0_slice: &[f64] = &component.k0.to_vec();
                    (material.re, material.im) =
                        regrid_lnk_data(l0_slice, n0_slice, k0_slice, &self.lam, true);
                    material.rho = component.rho;
                }
                RefractiveIndexKind::Other => {
                    let component = get_lnk_data(&material.key);
                    (material.re, material.im) =
                        regrid_lnk_data(component.l0, component.n0, component.k0, &self.lam, true);
                    material.rho = component.rho;
                }
            }
        }

        // Normalize the mass fractions
        (self.total_mfrac, self.tot_mfrac_core, self.tot_mfrac_mantle) = self
            .materials
            .iter()
            .fold((0.0, 0.0, 0.0), |(tot, tot_core, tot_mantle), m| {
                match m.kind {
                    MaterialKind::Core => (tot + m.mfrac, tot_core + m.mfrac, tot_mantle),
                    MaterialKind::Mantle => (tot + m.mfrac, tot_core, tot_mantle + m.mfrac),
                }
            });

        self.tot_mfrac_core /= self.total_mfrac;
        self.tot_mfrac_mantle /= self.total_mfrac;
        self.materials.iter_mut().for_each(|m| {
            m.mfrac /= self.total_mfrac;
        });
        self.materials.iter_mut().for_each(|m| {
            m.vfrac = m.mfrac / m.rho;
        });
        (self.tot_vfrac_core, self.tot_vfrac_mantle) =
            self.materials
                .iter()
                .fold((0.0, 0.0), |(tot_core, tot_mantle), m| match m.kind {
                    MaterialKind::Core => (tot_core + m.vfrac, tot_mantle),
                    MaterialKind::Mantle => (tot_core, tot_mantle + m.vfrac),
                });

        // Turn mass fractions into volume fractions
        self.rho_core = self.tot_mfrac_core / self.tot_vfrac_core;
        self.rho_mantle = self.tot_mfrac_mantle / self.tot_vfrac_mantle;

        // Normalize volume fractions
        self.materials.iter_mut().for_each(|m| match m.kind {
            MaterialKind::Core => m.vfrac /= self.tot_vfrac_core,
            MaterialKind::Mantle => m.vfrac /= self.tot_vfrac_mantle,
        });

        if self.pcore > 0.0 {
            self.materials.iter_mut().for_each(|m| {
                if m.kind == MaterialKind::Core {
                    m.vfrac *= 1.0 - self.pcore;
                }
            });
            self.rho_core *= 1.0 - self.pcore;
        }

        if self.pmantle > 0.0 {
            self.materials.iter_mut().for_each(|m| {
                if m.kind == MaterialKind::Mantle {
                    m.vfrac *= 1.0 - self.pmantle;
                }
            });
            self.rho_mantle *= 1.0 - self.pmantle;
        }

        // Calculate average density of the whole grain
        if self.nmant() == 0 {
            self.rho_av = self.rho_core;
        } else {
            self.rho_av = self.rho_core
                / (1.0 + self.tot_mfrac_mantle * (self.rho_core / self.rho_mantle - 1.0));
            self.tot_vfrac_mantle = self.tot_mfrac_mantle * self.rho_av / self.rho_mantle;
        }

        Ok(())
    }

    fn validate_size_distribution(&self) -> Result<(), KappaError> {
        match self.sizedis {
            SizeDistribution::PowerLaw if self.apow < 0.0 => Err(anyhow!("UnexpectedApow").into()),
            SizeDistribution::Normal => {
                if self.amean <= 0.0 {
                    return Err(
                        anyhow!("`amean` must be positive for (log-)normal distribution").into(),
                    );
                }
                if self.asigma == 0.0 {
                    return Err(
                        anyhow!("`asigma` cannot be zero for (log-)normal distribution").into(),
                    );
                }
                if self.asigma > 0.0 {
                    return Err(KappaError::ForceLogNormal);
                }
                Ok(())
            }
            _ => Ok(()),
        }
    }

    pub fn prepare_sparse(&mut self) {
        let f = self.lam[self.nlam - 1] / self.lam[self.nlam - 2];

        for (valmin, valmax) in self.scatlammin.iter_mut().zip(self.scatlammax.iter_mut()) {
            if *valmax / *valmin < f {
                *valmin /= f.sqrt() / 1.0001;
                *valmax *= f.sqrt() / 1.0001;
            }
        }

        for (il, &lam_val) in self.lam.iter().enumerate() {
            if self
                .scatlammin
                .iter()
                .zip(self.scatlammax.iter())
                .any(|(&min_val, &max_val)| lam_val >= min_val && lam_val <= max_val)
            {
                self.sparse_indices.insert(il);
                if il > 0 {
                    self.sparse_indices.insert(il - 1);
                }
                if il + 1 < self.nlam {
                    self.sparse_indices.insert(il + 1);
                }
            }
        }
    }

    /// Compute the moments of the size distribution
    ///
    /// The results are returned in `ameans`, an array of length 3:
    /// $$
    ///                [\langle a \rangle,\quad \langle a^2 \rangle^{1/2},\quad \langle a^3 \rangle^{1/3}]
    /// $$
    /// If both mn and sig are nonzero and the product is positive, we use
    /// the log-normal size distribution.  If not, we use the powerlaw.
    #[must_use]
    pub fn get_sizedis_moments(&self) -> [f64; 3] {
        let ns = 1000;

        let aminlog = self.amin.log10();
        let amaxlog = self.amax.log10();
        let pow = -self.apow;

        if ((self.amax - self.amin) / self.amin).abs() < 1e-6 {
            [self.amin, self.amin, self.amin]
        } else {
            let mut tot = [0.0; 3];
            let mut totn = 0.0;
            let mut ameans = [0.0; 3];
            for is in 0..ns {
                let isf = f64::from(is);
                let nsf = f64::from(ns);
                let r = 10.0f64.powf(aminlog + (amaxlog - aminlog) * isf / nsf);
                let nr = if (self.amean * self.asigma).abs() > 0.0 {
                    // normal or log-normal size distribution
                    let expo = if self.asigma > 0.0 {
                        0.5 * ((r - self.amean) / self.asigma).powi(2)
                    } else {
                        0.5 * ((r / self.amean).ln() / self.asigma).powi(2)
                    };
                    if expo > 99.0 { 0.0 } else { -expo.exp() }
                } else {
                    r.powf(pow + 1.0)
                };
                totn += nr;
                tot[0] += nr * r;
                tot[1] += nr * r.powi(2);
                tot[2] += nr * r.powi(3);
            }
            if totn == 0.0 {
                totn = 1.0;
            }
            ameans[0] = tot[0] / totn;
            ameans[1] = (tot[1] / totn).sqrt();
            ameans[2] = (tot[2] / totn).powf(1.0 / 3.0);
            ameans
        }
    }

    /// Write the wavelength grid to a file in the specified format.
    ///
    /// # Panics
    /// - If the file cannot be created or written to.
    /// - If the wavelength values are not positive (to accommodate log-log interpolation).
    /// # Errors
    /// - If the file format is invalid, such as missing header or malformed data lines.
    pub fn write_wavelength_grid(&self) -> Result<()> {
        let lamoutfile = self.outdir.clone() + "/kappa_lam.dat";
        let file = File::create(lamoutfile).expect("Failed to open file");
        let mut writer = BufWriter::new(file);

        writeln!(
            writer,
            "# Wavelength grid written by kappa, can be read back in with -l kappa_lam.dat"
        )?;
        writeln!(writer, "# First line: number of wavelengths")?;
        writeln!(writer, "# Then one lambda per line, in micrometer")?;
        writeln!(writer, "{}", self.nlam)?;
        for i in 0..self.nlam {
            writeln!(writer, "{:18.5e}", self.lam[i])?;
        }
        Ok(())
    }

    /// Run the simulation.
    ///
    /// # Panics
    /// - If the number of materials exceeds 20.
    /// - If the porosity values are not in the range [0, 1).
    /// - If the grain size distribution parameters are invalid.
    /// - If the wavelength grid parameters are invalid.
    /// - If the DHS `fmax` parameter is not in the range [0, 1).
    /// # Errors
    /// - `TooManyMaterials`: The number of materials exceeds 20.
    /// - `InvalidPorosity`: The porosity values are not in the range [0, 1).
    /// - `InvalidSizeInput`: The grain size distribution parameters are invalid.
    /// - `InvalidWavelengthInput`: The wavelength grid parameters are invalid.
    /// - `SamplingRequired`: More than one wavelength is required when `lmin` and `lmax` are different.
    pub fn run(&self) -> Result<(), KappaError> {
        // Loop for splitting the output into files by grain size
        if self.split {
            let mut nsub = self.nsubgrains;
            if nsub.is_multiple_of(2) {
                nsub += 1;
            }
            let afact = (self.amax / self.amin).powf(1.0 / self.na as f64);
            let afsub = afact.powf(1.0 / (nsub - 1) as f64);

            let bar = Arc::new(ProgressBar::new(self.nlam as u64));

            bar.set_style(
                ProgressStyle::default_bar()
                    .template(
                        "{spinner} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
                    )
                    .unwrap()
                    .progress_chars("█▓▒░ "),
            );

            (0..self.na).into_par_iter().try_for_each(|ia| {
                let iaf = ia as f64;
                let asplit = self.amin * afact.powf(iaf + 0.5);
                let nsubf = nsub as f64;
                let aminsplit = asplit * afsub.powf(-nsubf / 2.0);
                let amaxsplit = asplit * afsub.powf(nsubf / 2.0);
                let _ = self.compute_kappa(ia, aminsplit, amaxsplit)?;
                bar.inc(1); // Increment the progress bar inside the parallel loop
                Ok::<(), KappaError>(())
            })?;
            bar.finish_with_message("Processing complete!");
            exit(0);
        } else {
            let p = self.compute_kappa(0, self.amin, self.amax)?;
            if !p.is_ok {
                if self.mmf_ss {
                    eprintln!("WARNING: opacities OK, but some F_nn,g_asym may not be accurate");
                } else {
                    eprintln!("WARNING: opacities OK, but some F_nn,g_asym are set to zero");
                }
            }

            // TODO: FITS writer will be implemented later

            self.write_opacities("kappa_opacity", &p)?;
        }

        Ok(())
    }

    /// Calculate the opacities for a grain-size value.
    ///
    /// # Errors
    /// Returns a [`KappaError`] if there is an issue with parsing the arguments,
    /// such as missing required arguments, invalid argument formats, or unknown options.
    ///
    /// # Panics
    /// May panic if the number of materials exceeds 20, the porosity values are not in the
    /// appropriate range, or the grain size distribution parameters are invalid.
    pub fn compute_kappa(&self, ia: usize, amin: f64, amax: f64) -> Result<Particle, KappaError> {
        let ns = self.na;
        let (nf, ifmn) = if self.fmax == 0.0 { (1, 1) } else { (20, 12) };

        let mut r = RVector::zeros(ns);
        let mut nr = RVector::zeros(ns);

        let mut e1mantle = RVector::zeros(self.nlam);
        let mut e2mantle = RVector::zeros(self.nlam);

        let aminlog = amin.log10();
        let amaxlog = amax.log10();
        let pow = -self.apow;

        let mut tot = 0.0;

        if ns == 1 {
            // Just one size
            r[0] = 10.0_f64.powf(f64::midpoint(aminlog, amaxlog));
            nr[0] = r[0].powf(pow + 1.0); // should be 1/r[0]^3 ???  Not important.
        } else {
            // Size distribution
            let nsf = (ns - 1) as f64;
            for (is, (r_val, nr_val)) in r.iter_mut().zip(nr.iter_mut()).enumerate() {
                let isf = is as f64;
                *r_val = 10.0f64.powf(aminlog + (amaxlog - aminlog) * isf / nsf);
                // Apply size distribution and exponential calculation
                *nr_val = match self.sizedis {
                    // Power law size distribution
                    SizeDistribution::Normal | SizeDistribution::LogNormal => {
                        if self.asigma < 0.0 {
                            let expo = 0.5 * ((*r_val - self.amean) / self.asigma).powf(2.0);
                            if expo > 99.0 { 0.0 } else { (-expo).exp() }
                        } else {
                            let expo = 0.5 * ((*r_val / self.amean).ln() / self.asigma).powf(2.0);
                            if expo > 99.0 { 0.0 } else { (-expo).exp() }
                        }
                    }
                    SizeDistribution::File => {
                        todo!() // TODO: Move this out of else block
                    }
                    SizeDistribution::PowerLaw => r_val.powf(pow + 1.0),
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

        if self.write_grid && ia == 0 {
            self.write_sizedis_file(ns, &r, &mut nr, tot)?;
        }

        // Copy the refractory index data for all materials into local arrays
        let mut e1 = RMatrix::zeros((self.nlam, self.materials.len()));
        let mut e2 = RMatrix::zeros((self.nlam, self.materials.len()));
        for im in 0..self.materials.len() {
            e1.slice_mut(s![.., im]).assign(&self.materials[im].re);
            e2.slice_mut(s![.., im]).assign(&self.materials[im].im);
        }

        let mut e1_blend = RVector::zeros(self.nlam);
        let mut e2_blend = RVector::zeros(self.nlam);
        let mut e_in = CVector::zeros(self.ncore());

        // Mixing, for all wavelengths
        for il in 0..self.nlam {
            // Core
            if self.materials.len() == 1 && self.pcore == 0.0 {
                // Solid core, single material, nothing to blend for the core
                e1_blend[il] = e1[[il, 0]];
                e2_blend[il] = e2[[il, 0]];
            } else {
                // Blend the core materials
                for im in 0..self.ncore() {
                    e_in[im] = Complex::new(e1[[il, im]], e2[[il, im]]);
                }
                let e_out = bruggeman_blend(
                    &self
                        .materials
                        .iter()
                        .filter(|m| matches!(m.kind, MaterialKind::Core))
                        .map(|m| m.vfrac)
                        .collect::<Vec<_>>(),
                    &e_in.to_vec(),
                )?;
                e1_blend[il] = e_out.re;
                e2_blend[il] = e_out.im;
            }
            // dbg!(&e1_blend[il], &e2_blend[il]);
            // Mantle
            if self.nmant() > 0 {
                // We do have a mantle to add
                if (self.nmant() == 1) && (self.pmantle == 0.0) {
                    // No Blending needed inside the mantle - just copy e1 and e2
                    // Since it is only one material, we know it is index nm
                    e1mantle[il] = e1[[il, self.materials.len() - 1]];
                    e2mantle[il] = e2[[il, self.materials.len() - 1]];
                } else {
                    // Blend the mantle materials
                    for im in 0..self.nmant() {
                        e_in[im] =
                            Complex::new(e1[[il, im + self.ncore()]], e2[[il, im + self.ncore()]]);
                    }
                    let e_out = bruggeman_blend(
                        &self
                            .materials
                            .iter()
                            .filter(|m| matches!(m.kind, MaterialKind::Mantle))
                            .map(|m| m.vfrac)
                            .collect::<Vec<_>>(),
                        &e_in.to_vec(),
                    )?;
                    e1mantle[il] = e_out.re;
                    e2mantle[il] = e_out.im;
                }
                let (e1mg, e2mg) = maxwell_garnet_blend(
                    Complex::new(e1_blend[il], e2_blend[il]), // Core
                    Complex::new(e1mantle[il], e2mantle[il]), // Mantle
                    self.tot_vfrac_mantle,
                );
                e1_blend[il] = e1mg;
                e2_blend[il] = e2mg;
            }
        }
        // dbg!(&e1_blend, &e2_blend);
        // exit(1);

        if self.blend_only {
            // Write .lnk file
            todo!()
        }
        // TODO: Now there is some `--print` command logic
        // TODO: that we will implement later

        if self.write_grid || self.blend_only {
            println!("Exiting after writing requested files");
            exit(0);
        }

        // Check how we are going to average over hollow sphere components
        let (f, wf) = if nf > 1 && self.fmax > 0.01 {
            // Get the weights for Gauss-Legendre integration
            gauss_legendre(0.01f64, self.fmax, nf)
        } else if self.fmax == 0.0 {
            (RVector::zeros(nf), RVector::ones(nf))
        } else {
            // Just a compact sphere, weight is 1
            (RVector::ones(nf) * self.fmax, RVector::ones(nf))
        };

        // Initialize mu
        let mut mu = RVector::zeros(self.nang);
        let nangby2f = (self.nang / 2) as f64;
        for (j, val) in mu.iter_mut().take(self.nang / 2).enumerate() {
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

        let p =
            {
                let bar = ProgressBar::new(self.nlam as u64);
                bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap()
                .progress_chars("█▓▒░ "),
        );
                let counter = Arc::new(AtomicUsize::new(0));
                let results: Result<Vec<KappaResult>, _> = if self.split {
                    (0..self.nlam)
                        .map(|ilam| self.over_wavelengths(ilam, &kps))
                        .collect()
                } else {
                    let counter = Arc::clone(&counter);
                    let bar = bar.clone();
                    let res = (0..self.nlam)
                        .into_par_iter()
                        .map(|ilam| {
                            let r = self.over_wavelengths(ilam, &kps);
                            let prev = counter.fetch_add(1, AtomicOrdering::SeqCst);
                            bar.set_position((prev + 1) as u64);
                            r
                        })
                        .collect();
                    bar.finish_with_message("Processing complete!");
                    res
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

                    for ilam in 0..self.nlam {
                        if !trust[ilam] {
                            is_ok = false;
                            is_ok_lmin = self.lam[ilam];
                        }
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

    /// Calculate the opacities for a wavelength value.
    fn over_wavelengths(&self, ilam: usize, kps: &KappaState) -> Result<KappaResult> {
        let qabsdqext_min = 1e-4;
        let wvno = 2.0 * PI / self.lam[ilam];
        let mut f11 = RVector::zeros(self.nang);
        let mut f12 = RVector::zeros(self.nang);
        let mut f22 = RVector::zeros(self.nang);
        let mut f33 = RVector::zeros(self.nang);
        let mut f34 = RVector::zeros(self.nang);
        let mut f44 = RVector::zeros(self.nang);
        let mut csca = 0.0;
        let mut cabs = 0.0;
        let mut cext = 0.0;
        let mut smat_nbad = 0;
        let mut mass = 0.0;
        let mut vol = 0.0;
        let nangf = self.nang as f64;
        let mut mie_f11 = RVector::zeros(self.nang);
        let mut mie_f12 = RVector::zeros(self.nang);
        let mut mie_f22 = RVector::zeros(self.nang);
        let mut mie_f33 = RVector::zeros(self.nang);
        let mut mie_f34 = RVector::zeros(self.nang);
        let mut mie_f44 = RVector::zeros(self.nang);
        for is in 0..kps.ns {
            let r1 = kps.r[is];
            match self.method {
                KappaMethod::DHS => {
                    let m_in = Complex::new(1.0, 0.0);
                    for ifn in 0..kps.nf {
                        let mut rad: f64;
                        let mut spheres = false;
                        let mut too_large = false;
                        let mut scat_mie: f64 = 0.0;
                        let mut ext_mie: f64 = 0.0;
                        if kps.f[ifn] == 0.0 {
                            spheres = true;
                        } else if r1 * wvno > self.xlim {
                            too_large = true;
                        } else {
                            rad = r1 / (1.0 - kps.f[ifn]).powf(1.0 / 3.0);
                            let rcore = rad * kps.f[ifn].powf(1.0 / 3.0);
                            let mconj = Complex::new(kps.e1_blend[ilam], -kps.e2_blend[ilam]);
                            let wvno_1 = wvno.min(self.xlim_dhs / rad);
                            let dhs_cfg = DHSConfig {
                                r_core: rcore,
                                r_shell: rad,
                                wave_number: wvno_1,
                                r_indsh: mconj,
                                r_indco: m_in,
                                mu: kps.mu,
                                numang: self.nang / 2,
                                max_angle: self.nang,
                            };

                            match toon_ackerman_1981(&dhs_cfg) {
                                Ok(dhs_result) => {
                                    // DHS succeeded
                                    let ext_mie_factor = dhs_result.q_ext * PI * rad.powi(2);
                                    let scat_mie_factor = dhs_result.q_sca * PI * rad.powi(2);
                                    let factor = 2.0 * PI / scat_mie_factor / wvno.powi(2);
                                    for j in 0..(self.nang / 2) {
                                        mie_f11[j] = (dhs_result.m1[[j, 0]]
                                            + dhs_result.m0[[j, 0]])
                                            * factor;
                                        mie_f12[j] = (dhs_result.m1[[j, 0]]
                                            - dhs_result.m0[[j, 0]])
                                            * factor;
                                        mie_f22[j] = (dhs_result.m1[[j, 0]]
                                            + dhs_result.m0[[j, 0]])
                                            * factor;
                                        mie_f33[j] = dhs_result.s10[[j, 0]] * factor;
                                        mie_f34[j] = -dhs_result.d10[[j, 0]] * factor;
                                        mie_f44[j] = dhs_result.s10[[j, 0]] * factor;
                                        mie_f11[self.nang - j - 1] = (dhs_result.m1[[j, 1]]
                                            + dhs_result.m0[[j, 1]])
                                            * factor;
                                        mie_f12[self.nang - j - 1] = (dhs_result.m1[[j, 1]]
                                            - dhs_result.m0[[j, 1]])
                                            * factor;
                                        mie_f22[self.nang - j - 1] = (dhs_result.m1[[j, 1]]
                                            + dhs_result.m0[[j, 1]])
                                            * factor;
                                        mie_f33[self.nang - j - 1] =
                                            dhs_result.s10[[j, 1]] * factor;
                                        mie_f34[self.nang - j - 1] =
                                            -dhs_result.d10[[j, 1]] * factor;
                                        mie_f44[self.nang - j - 1] =
                                            dhs_result.s10[[j, 1]] * factor;
                                    }
                                    scat_mie = scat_mie_factor;
                                    ext_mie = ext_mie_factor;
                                }
                                Err(e) => {
                                    println!("DHS error: {e:?}");
                                    // err==1: use compact sphere rad=r1
                                    // rad = r1;
                                    spheres = true; // reuse the Mie path below
                                }
                            }
                        }

                        // Mie fallback for: spheres, too_large, or DHS error (spheres flag set on err)
                        if spheres || too_large {
                            if too_large {
                                rad = r1 / (1.0 - kps.f[kps.ifmn]).powf(1.0 / 3.0);
                            } else {
                                rad = r1;
                            }
                            let lmie = self.lam[ilam];
                            let e1_mie = kps.e1_blend[ilam];
                            let e2_mie = kps.e2_blend[ilam];
                            let thmin = 180.0 * 0.5 / nangf;
                            let thmax = 180.0 * (nangf - 0.5) / nangf;
                            let radius = if rad / lmie < 5000.0 {
                                rad
                            } else {
                                rad / 5000.0
                            };
                            let miec = MieConfig {
                                nangle: self.nang,
                                delta: 1e-8,
                                thmin,
                                thmax,
                                step: (thmax - thmin) / (nangf - 1.0),
                                lam: lmie,
                                cmm: Complex::new(e1_mie, e2_mie),
                                rad: radius,
                            };
                            let mie_result = de_rooij_1984(&miec)?;
                            scat_mie = mie_result.c_sca;
                            ext_mie = mie_result.c_ext;
                            // populate matrix from mie_result (adjust field names as needed)
                            for j in 0..self.nang {
                                mie_f11[j] = mie_result.f_11[j];
                                mie_f12[j] = mie_result.f_12[j];
                                mie_f22[j] = mie_result.f_11[j];
                                mie_f33[j] = mie_result.f_33[j];
                                mie_f34[j] = mie_result.f_34[j];
                                mie_f44[j] = mie_result.f_33[j];
                            }
                        }

                        let mut tot = 0.0;
                        let mut tot2 = 0.0;
                        for j in 0..self.nang {
                            let jf = j as f64;
                            tot += mie_f11[j] * (PI * (jf + 0.5) / nangf).sin();
                            tot2 += (PI * (jf + 0.5) / nangf).sin();
                        }
                        mie_f11[0] += (tot2 - tot) / (PI * 0.5 / nangf).sin();
                        if mie_f11[0] < 0.0 {
                            mie_f11[0] = 0.0;
                        }
                        for j in 0..self.nang {
                            f11[j] += kps.wf[ifn] * kps.nr[is] * mie_f11[j] * scat_mie;
                            f12[j] += kps.wf[ifn] * kps.nr[is] * mie_f12[j] * scat_mie;
                            f22[j] += kps.wf[ifn] * kps.nr[is] * mie_f22[j] * scat_mie;
                            f33[j] += kps.wf[ifn] * kps.nr[is] * mie_f33[j] * scat_mie;
                            f34[j] += kps.wf[ifn] * kps.nr[is] * mie_f34[j] * scat_mie;
                            f44[j] += kps.wf[ifn] * kps.nr[is] * mie_f44[j] * scat_mie;
                        }
                        let mut camie = ext_mie - scat_mie;
                        if camie < ext_mie * qabsdqext_min {
                            camie = ext_mie * 1e-4;
                            ext_mie = camie + scat_mie;
                        }
                        cext += kps.wf[ifn] * kps.nr[is] * ext_mie;
                        csca += kps.wf[ifn] * kps.nr[is] * scat_mie;
                        cabs += kps.wf[ifn] * kps.nr[is] * camie;
                        mass +=
                            kps.wf[ifn] * kps.nr[is] * self.rho_av * 4.0 * PI * r1.powi(3) / 3.0;
                        vol += kps.wf[ifn] * kps.nr[is] * 4.0 * PI * r1.powi(3) / 3.0;
                    }
                }
                KappaMethod::MMF => {
                    let m_mono = 4.0 * PI / 3.0 * self.mmf_a0.powi(3) * self.rho_av;
                    let v_agg = 4.0 * PI / 3.0 * r1.powi(3);
                    let m_agg = v_agg * self.rho_av;
                    let nmono = m_agg / m_mono;
                    let dfrac = if self.mmf_struct > 1.0 {
                        self.mmf_struct
                    } else {
                        3.0 * nmono.ln() / (nmono - self.mmf_struct).ln()
                    };
                    let kfrac = if self.mmf_kf > 0.0 {
                        self.mmf_kf
                    } else {
                        (5.0 / 3.0).powf(dfrac)
                    };
                    let fracc = FractalConfig {
                        solver: FractalSolver::ModifiedMeanField,
                        cutoff: FractalCutoff::Gaussian,
                        geometry: FractalGeometry::Tazaki,
                        nang: self.nang / 2 + 1,
                        pn: nmono,
                        r0: self.mmf_a0,
                        df: dfrac,
                        k0: kfrac,
                        lmd: self.lam[ilam],
                        refrel: Complex::new(kps.e1_blend[ilam], kps.e2_blend[ilam]),
                        ..Default::default()
                    };
                    let frac_result = mean_scattering(&fracc)?;
                    if frac_result.dphi > 1.0 {
                        smat_nbad += 1;
                    }
                    let factor = 4.0 * PI / wvno.powi(2);
                    for j in 0..self.nang {
                        f11[j] += kps.nr[is]
                            * 0.5
                            * (frac_result.smat[[0, j]] + frac_result.smat[[0, j + 1]])
                            * factor;
                        f12[j] += kps.nr[is]
                            * 0.5
                            * (frac_result.smat[[1, j]] + frac_result.smat[[1, j + 1]])
                            * factor;
                        f22[j] += kps.nr[is]
                            * 0.5
                            * (frac_result.smat[[0, j]] + frac_result.smat[[0, j + 1]])
                            * factor;
                        f33[j] += kps.nr[is]
                            * 0.5
                            * (frac_result.smat[[2, j]] + frac_result.smat[[2, j + 1]])
                            * factor;
                        f34[j] += kps.nr[is]
                            * 0.5
                            * (frac_result.smat[[3, j]] + frac_result.smat[[3, j + 1]])
                            * factor;
                        f44[j] += kps.nr[is]
                            * 0.5
                            * (frac_result.smat[[2, j]] + frac_result.smat[[2, j + 1]])
                            * factor;
                    }
                    cext += kps.nr[is] * frac_result.c_ext;
                    csca += kps.nr[is] * frac_result.c_sca;
                    cabs += kps.nr[is] * frac_result.c_abs;
                    mass += kps.nr[is] * m_agg;
                    vol += kps.nr[is] * v_agg;
                }
                KappaMethod::CDE => {
                    let v = 4.0 * PI * r1.powi(3) / 3.0;
                    vol += kps.nr[is] * v;
                    mass += kps.nr[is] * self.rho_av * v;
                    let m = Complex::new(kps.e1_blend[ilam], kps.e2_blend[ilam]);
                    let dcabs =
                        2.0 * wvno * v * (m.powi(2) / (m.powi(2) - 1.0) * m.powi(2).ln()).im;
                    let dcsca = if (kps.e2_blend[ilam] / kps.e1_blend[ilam]).abs() < 1e-6 {
                        // Non-absorbing case
                        let mre = m.re;
                        wvno.powi(4) * v.powi(2) / (3.0 * PI)
                            * (mre.powi(2) - 1.0 - mre.powi(2).ln())
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
                        let lmie = self.lam[ilam];
                        let rmie = lmie / 1e3;
                        let e1_mie = kps.e1_blend[self.nlam - 1];
                        let e2_mie = kps.e2_blend[self.nlam - 1];
                        let thmin = 180.0 * (1.0 - 0.5) / nangf;
                        let thmax = 180.0 * (nangf - 0.5) / nangf;
                        let miec = MieConfig {
                            nangle: self.nang,
                            delta: 1e-8,
                            thmin,
                            thmax,
                            step: (thmax - thmin) / (nangf - 1.0),
                            lam: lmie,
                            cmm: Complex::new(e1_mie, e2_mie),
                            rad: rmie,
                        };
                        let mie_result = de_rooij_1984(&miec)?;
                        let mie_f22 = mie_result.f_11.clone();
                        let mie_f44 = mie_result.f_33.clone();
                        for j in 0..self.nang {
                            f11[j] = mie_result.f_11[j] * csca;
                            f12[j] = mie_result.f_12[j] * csca;
                            f22[j] = mie_f22[j] * csca;
                            f33[j] = mie_result.f_33[j] * csca;
                            f34[j] = mie_result.f_34[j] * csca;
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

        if self.chop_angle > 0.0 {
            // Flatten the peak
            let ichop = (self.chop_angle / (180.0 / nangf)) as usize;
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
            for j in 0..self.nang {
                let jf = j as f64;
                tot += f11[j] * (PI * (jf + 0.5) / nangf).sin() * (PI / nangf) * 2.0 * PI;
                tot2 += (PI * (jf + 0.5) / nangf).sin() * (PI / nangf) * 2.0 * PI;
            }
            for j in 0..self.nang {
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
        for i in 0..self.nang {
            let ir = i as f64;
            g += f11[i] * (PI * (ir + 0.5) / nangf).cos() * (PI * (ir + 0.5) / nangf).sin();
            tot += f11[i] * (PI * (ir + 0.5) / nangf).sin();
        }
        g /= tot;
        if smat_nbad > 0 && !self.mmf_ss {
            g = 0.0;
            for i in 0..self.nang {
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

    fn create_output_file(&self, filename: &str) -> Result<File> {
        // Create directory if it doesn't exist
        fs::create_dir_all(&self.outdir)?;

        // Construct full file path
        let filepath = format!("{}/{}", self.outdir, filename);

        // Remove existing file if it exists
        if Path::new(&filepath).exists() {
            fs::remove_file(&filepath)?;
        }

        // Create and return the new file handle
        File::create(&filepath).map_err(|_| anyhow!("Failed to create file!"))
    }

    /// Write the computed opacities to a file in the specified format.
    ///
    /// # Panics
    /// - If the file cannot be created or written to.
    /// - If the size distribution values are not positive (to accommodate log-log interpolation).
    /// # Errors
    /// - If the file format is invalid, such as missing header or malformed data lines.
    pub fn write_opacities(&self, fname: &str, p: &Particle) -> Result<()> {
        let (ext, _ml) = if self.for_radmc {
            unimplemented!()
        } else {
            ("dat", "F11 F12 F22 F33 F34 F44")
        };

        let opacity_filename = format!("{fname}.{ext}");

        let file_opacity = self.create_output_file(&opacity_filename)?;
        let mut writer = BufWriter::new(file_opacity);

        // Write header
        let cc = "#";
        let ameans = self.get_sizedis_moments();
        writeln!(
            writer,
            "{cc}============================================================================"
        )?;
        let a_vals = match self.sizedis {
            SizeDistribution::File => self.ameans_file,
            _ => ameans,
        };
        writeln!(
            writer,
            "{} Opacities computed by kappa        <a^n>= {:11.4e} {:11.4e} {:11.4e}",
            cc, a_vals[0], a_vals[1], a_vals[2]
        )?;
        match self.method {
            KappaMethod::MMF => {
                let s_struct = if self.mmf_struct > 1.0 {
                    "(fractal dimension)"
                } else {
                    "(filling factor)"
                };
                writeln!(
                    writer,
                    "{} Method:   {:?}  a0={:7.3}  Struct={:7.3}{}",
                    cc, self.method, self.mmf_a0, self.mmf_struct, s_struct
                )?;
            }
            KappaMethod::CDE => {
                writeln!(writer, "{} Method:   {:?}", cc, self.method)?;
            }
            KappaMethod::DHS => {
                writeln!(
                    writer,
                    "{} Method:   {:?}  fmax={:7.3}",
                    cc, self.method, self.fmax
                )?;
            }
        }
        writeln!(writer, "{cc} Parameters:")?;
        match self.sizedis {
            SizeDistribution::File => {
                // writeln!(
                //     writer,
                //     "{}   amin [um]={:11.3} amax [um]={:11.3}  na  ={:5}    file={}",
                //     cc,
                //     kpc.amin,
                //     kpc.amax,
                //     kpc.na,
                //     kpc.
                // )?;
                todo!()
            }
            SizeDistribution::LogNormal | SizeDistribution::Normal => {
                writeln!(
                    writer,
                    "{}   amin [um]={:11.3} amax [um]={:11.3}  na  ={:5}    {:?}={:.4}:{:.4}",
                    cc, self.amin, self.amax, self.na, self.sizedis, self.amean, self.asigma
                )?;
            }
            SizeDistribution::PowerLaw => {
                writeln!(
                    writer,
                    "{}   amin [um]={:11.3} amax [um]={:11.3}  na  ={:5}     apow={:10.2}",
                    cc, self.amin, self.amax, self.na, self.apow
                )?;
            }
        }
        writeln!(
            writer,
            "{}   lmin [um]={:11.3} lmax [um]={:11.3}  nlam={:5}     nang={:6}",
            cc, self.lmin, self.lmax, self.nlam, self.nang
        )?;
        writeln!(
            writer,
            "{}   porosity ={:11.3} p_mantle ={:11.3}  fmax={:9.2} chop=  {:4.1}",
            cc, self.pcore, self.pmantle, self.fmax, self.chop_angle
        )?;

        writeln!(writer, "{cc} Composition:")?;
        writeln!(writer, "{cc}  Where   mfrac  rho   Material")?;
        writeln!(
            writer,
            "{cc}  -----   -----  ----  -----------------------------------------------------"
        )?;

        let total_mfrac: f64 = self.materials.iter().map(|m| m.mfrac).sum();
        for mat in &self.materials {
            writeln!(
                writer,
                "{}  {:?} {:7.3} {:6.2}  {}",
                cc,
                mat.kind,
                mat.mfrac / total_mfrac,
                mat.rho,
                mat.key,
            )?;
        }

        if self.rho_av > 0.0 {
            writeln!(
                writer,
                "{cc}  - - -   - - -  -  -  - - - - - - - - - - - - - - - - - - - - - - - - - - -",
            )?;
            let label = if (self.pcore + self.pmantle) > 0.0 {
                "materials and vacuum"
            } else {
                "materials"
            };
            writeln!(
                writer,
                "{}  {:<6} {:7.3} {:6.2}  mixture of {:3} {}",
                cc,
                "grain",
                1.0,
                self.rho_av,
                self.materials.len(),
                label
            )?;
        }

        writeln!(
            writer,
            "{cc}----------------------------------------------------------------------------"
        )?;
        // writeln!(writer, "{} Command: {}", cc, kpc.)?;
        writeln!(
            writer,
            "{cc}============================================================================"
        )?;
        // end of header

        // write the output
        let header_line = if self.for_radmc {
            "# Output file formatted for RADMC-3D, dustkappa, no scattering matrix"
        } else {
            "# Standard output file, no scattering matrix"
        };
        writeln!(writer, "{header_line}")?;
        writeln!(writer, "#    iformat")?;
        writeln!(writer, "#    nlambda")?;
        writeln!(
            writer,
            "#    lambda[um]  kabs [cm^2/g]  ksca [cm^2/g]    g_asymmetry"
        )?;
        writeln!(
            writer,
            "#============================================================================"
        )?;
        writeln!(writer, "{}", 3)?; // iformat
        writeln!(writer, "{}", self.nlam)?; // nlambda

        for i in 0..self.nlam {
            writeln!(
                writer,
                "{:15.6e} {:15.6e} {:15.6e} {:15.6e}",
                self.lam[i], p.k_abs[i], p.k_sca[i], p.g[i]
            )?;
        }

        writer.flush()?;

        Ok(())
    }

    /// Write the size distribution to a file in the specified format.
    ///
    /// # Panics
    /// - If the file cannot be created or written to.
    /// - If the size distribution values are not positive (to accommodate log-log interpolation).
    /// # Errors
    /// - If the file format is invalid, such as missing header or malformed data lines.
    pub fn write_sizedis_file(
        &self,
        ns: usize,
        r: &RVector,
        nr: &mut RVector,
        tot: f64,
    ) -> Result<()> {
        let sdoutfile = self.outdir.clone() + "/kappa_sd.dat";
        let file = File::create(sdoutfile).expect("Failed to open file");
        let mut writer = BufWriter::new(file);

        writeln!(
            writer,
            "# Size distribution written by kappa, can be read in with -a kappa_sd.dat"
        )?;
        if self.split {
            writeln!(
                writer,
                "# This is only the first subparticle because of the -d switch"
            )?;
        }
        writeln!(writer, "# First line: Number of grain size bins NA")?;
        writeln!(writer, "# Then NA lines with:  agrain[um]  n(a)")?;
        writeln!(writer, "#   n(a) is the number of grains in the bin.")?;
        writeln!(
            writer,
            "#   In a linear grid      (da   =const), this would be n(a) = f(a)*da"
        )?;
        writeln!(
            writer,
            "#   In a logarithmic grid (dloga=const), this would be n(a) = f(a)*a*dloga."
        )?;
        writeln!(
            writer,
            "#   In an arbitrary grid,  just give the number of grains in the bin."
        )?;
        writeln!(
            writer,
            "# No normalization is necessary, it is done automatically."
        )?;
        writeln!(writer, "{ns}")?;

        for i in 0..ns {
            nr[i] /= tot;
            writeln!(writer, "{:18.5e} {:18.5e}", r[i], nr[i])?;
        }
        Ok(())
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
        kpc.materials.push(Material {
            key: "h2o-w".into(),
            kind: MaterialKind::Core,
            mfrac: 0.2000,
            ..Default::default()
        });
        kpc
    }
}

#[derive(Debug)]
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

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize, Clone))]
#[derive(PartialEq, Debug)]
pub enum SizeDistribution {
    PowerLaw,
    File,
    Normal,
    LogNormal,
}

#[derive(Debug)]
pub enum WavelengthKind {
    Other,
    File,
}

impl SizeDistribution {}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize, Clone))]
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

/// Converts a [`crate::components::StaticComponent`] by value into a [`Component`].
///
/// This consumes the input [`crate::components::StaticComponent`], moving its data into the new [`Component`].
impl From<crate::components::StaticComponent> for Component {
    fn from(sc: crate::components::StaticComponent) -> Self {
        Component {
            name: sc.name.to_string(),
            class: sc.class.to_string(),
            state: sc.state.to_string(),
            rho: sc.rho,
            size: sc.size,
            l0: RVector::from(sc.l0.to_vec()),
            n0: RVector::from(sc.n0.to_vec()),
            k0: RVector::from(sc.k0.to_vec()),
        }
    }
}

/// Converts a reference to a [`crate::components::StaticComponent`] into a [`Component`].
///
/// This clones data from the borrowed [`crate::components::StaticComponent`], allowing reuse of the original.
impl From<&crate::components::StaticComponent> for Component {
    fn from(sc: &crate::components::StaticComponent) -> Self {
        Component {
            name: sc.name.to_string(),
            class: sc.class.to_string(),
            state: sc.state.to_string(),
            rho: sc.rho,
            size: sc.size,
            l0: RVector::from(sc.l0.to_vec()),
            n0: RVector::from(sc.n0.to_vec()),
            k0: RVector::from(sc.k0.to_vec()),
        }
    }
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
