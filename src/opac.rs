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
}
