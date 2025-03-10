//! The heart of `kappa`.
//!
//!

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
