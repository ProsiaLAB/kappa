//! Command line interface for the `kappa`.
use std::env;
use std::iter::Peekable;
use std::path::Path;

use anyhow::Result;
use colored::{Color, Colorize};

// use crate::io::read_lnk_file;
use crate::components::MATERIAL_KEYS;
use crate::io::{read_lnk_file, read_sizedis_file};
use crate::opac::{KappaConfig, KappaError, KappaMethod, SpecialConfigs};

enum SizeArg {
    PowerLaw {
        amin: f64,
        amax: Option<f64>,
        apow: Option<f64>,
        na: Option<usize>,
    },
    LogNormal {
        amin: f64,
        amax: f64,
        amean: f64,
        asigma: f64,
        na: Option<usize>,
    },
    File(String),
}

#[derive(Debug)]
pub enum MaterialType {
    Core,
    Mantle,
}

pub fn launch() -> Result<KappaConfig, KappaError> {
    let mut kpc = KappaConfig::default();

    let mut args = env::args().skip(1).peekable();
    if args.peek().is_none() {
        return Ok(kpc);
    }

    let mut materials: Vec<String> = Vec::with_capacity(20);
    let mut materials_mfrac: Vec<f64> = Vec::with_capacity(20);
    let mut material_types: Vec<MaterialType> = Vec::with_capacity(20);
    let mut materials_n: Vec<f64> = Vec::with_capacity(20);
    let mut materials_k: Vec<f64> = Vec::with_capacity(20);
    let mut materials_rho: Vec<f64> = Vec::with_capacity(20);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-c" => {
                if let Some(material) = args.next() {
                    if let Some((n, k, rho)) = process_material(
                        &material.to_string(),
                        &mut args,
                        &mut materials,
                        &mut materials_mfrac,
                        &mut material_types,
                        MaterialType::Core,
                    )? {
                        materials_n.push(n);
                        materials_k.push(k);
                        materials_rho.push(rho);
                    } else {
                        eprintln!("Missing material key/file path after -c/-m");
                        return Err(KappaError::MissingArgument("-c/-m".to_string()));
                    };
                } else {
                    eprintln!("Missing material key/file path after -c/-m");
                    return Err(KappaError::MissingArgument("-c/-m".to_string()));
                }
            }
            "-m" => {
                // Same as "-c"
                if let Some(material) = args.next() {
                    if let Some((n, k, rho)) = process_material(
                        &material.to_string(),
                        &mut args,
                        &mut materials,
                        &mut materials_mfrac,
                        &mut material_types,
                        MaterialType::Mantle,
                    )? {
                        materials_n.push(n);
                        materials_k.push(k);
                        materials_rho.push(rho);
                    } else {
                        eprintln!("Missing material key/file path after -c/-m");
                        return Err(KappaError::MissingArgument("-c/-m".to_string()));
                    };
                } else {
                    eprintln!("Missing material key/file path after -c/-m");
                    return Err(KappaError::MissingArgument("-c/-m".to_string()));
                }
            }
            "--diana" => {
                kpc = KappaConfig::diana();
            }
            "--dsharp" => {
                kpc = KappaConfig::dsharp();
            }
            "--dsharp-no-ice" => {
                kpc = KappaConfig::dsharp_no_ice();
            }
            // Size distribution options
            "-a" => match process_grain_size(&mut args)? {
                SizeArg::PowerLaw {
                    amin,
                    amax,
                    apow,
                    na,
                } => {
                    kpc.amin = amin;
                    if let Some(amax) = amax {
                        kpc.amax = amax;
                    }
                    if let Some(apow) = apow {
                        kpc.apow = apow;
                    }
                    if let Some(na) = na {
                        kpc.na = na;
                    }
                }
                SizeArg::LogNormal {
                    amin,
                    amax,
                    amean,
                    asigma,
                    na,
                } => {
                    kpc.amin = amin;
                    kpc.amax = amax;
                    kpc.amean = amean;
                    kpc.asigma = asigma;
                    if let Some(na) = na {
                        kpc.na = na;
                    }
                }
                SizeArg::File(file) => {
                    // Read file
                    let _ = read_sizedis_file(&file);
                }
            },
            "--amin" => {
                kpc.amin = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--amin".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--amax" => {
                kpc.amax = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--amax".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--amean" => {
                kpc.amean = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--amean".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--asig" | "--asigma" => {
                kpc.asigma = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--asig".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--apow" => {
                kpc.apow = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--apow".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--na" => {
                kpc.na = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--na".to_string()))?
                    .parse::<usize>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            // Wavelength options
            "-l" => {
                todo!()
            }
            "--lmin" => {
                kpc.lmin = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--lmin".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--lmax" => {
                kpc.lmax = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--lmax".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--nlam" => {
                kpc.nlam = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--nlam".to_string()))?
                    .parse::<usize>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            // Grain geometry and computational method options
            "-p" | "--porosity" => {
                kpc.pcore = args
                    .next()
                    .ok_or_else(|| KappaError::MissingArgument("--porosity".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
                kpc.pmantle = match args.next() {
                    Some(arg) => arg.parse::<f64>().map_err(|_| KappaError::InvalidType)?,
                    None => kpc.pcore, // Default is to use the same porosity
                };
            }
            "--dhs" | "--fmax" | "--mie" => {
                kpc.method = KappaMethod::DHS;
                if arg == "--mie" {
                    kpc.fmax = 0.0;
                } else if arg == "--fmax" {
                    // Peek at the next argument to check if it's a valid number
                    if let Some(next_arg) = args.peek() {
                        if let Ok(fmax_value) = next_arg.parse::<f64>() {
                            args.next(); // Consume the number
                            kpc.fmax = fmax_value;
                        } else {
                            kpc.fmax = 0.8; // Default value when no valid number follows
                        }
                    } else {
                        kpc.fmax = 0.8; // Default value when --fmax is the last argument
                    }
                } else {
                    kpc.fmax = 0.8; // Default value if neither --mie nor --fmax is provided
                }
            }
            "--xlim" => {
                kpc.xlim = args
                    .next()
                    .ok_or_else(|| KappaError::InvalidArgument("--xlim".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--xlim-dhs" => {
                kpc.xlim_dhs = args
                    .next()
                    .ok_or_else(|| KappaError::InvalidArgument("--xlim-dhs".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--mmf" | "--mmf-ss" => {
                kpc.method = KappaMethod::MMF;
                if arg == "--mmf-ss" {
                    kpc.mmf_ss = args
                        .next()
                        .ok_or_else(|| KappaError::InvalidArgument("--mmf-ss".to_string()))?
                        .parse::<bool>()
                        .map_err(|_| KappaError::InvalidType)?;
                }
                kpc.mmf_a0 = args
                    .next()
                    .ok_or_else(|| KappaError::InvalidArgument("--mmf-a0".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
                kpc.mmf_struct = args
                    .next()
                    .ok_or_else(|| KappaError::InvalidArgument("--mmf-struct".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
                kpc.mmf_kf = args
                    .next()
                    .ok_or_else(|| KappaError::InvalidArgument("--mmf-kf".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--cde" => {
                kpc.method = KappaMethod::CDE;
            }
            // Miscellaneous options
            "-L" | "--list" => {
                list_materials();
                return Ok(kpc);
            }
            "-h" => {
                print_short_help();
                return Ok(kpc);
            }
            "--help" => {
                print_long_help();
                return Ok(kpc);
            }
            "-q" | "--quiet" => {
                todo!()
            }
            "-v" | "--version" => {
                todo!()
            }
            // Unknown option
            _ if arg.starts_with('-') => {
                println!("Unknown option: {}", arg);
                return Ok(kpc);
            }
            // Positional arguments
            _ => {
                if let Some((n, k, rho)) = process_material(
                    arg.as_str(),
                    &mut args,
                    &mut materials,
                    &mut materials_mfrac,
                    &mut material_types,
                    MaterialType::Core,
                )? {
                    materials_n.push(n);
                    materials_k.push(k);
                    materials_rho.push(rho);
                } else {
                    eprintln!("Missing material key/file path after -c/-m");
                    return Err(KappaError::MissingArgument("-c/-m".to_string()));
                }
            }
        }
    }
    println!("args: {:?}", materials);
    println!("args: {:?}", materials_mfrac);
    println!("args: {:?}", material_types);
    // todo!()
    Ok(kpc)
}

fn process_grain_size<I>(args: &mut I) -> Result<SizeArg, KappaError>
where
    I: Iterator<Item = String>,
{
    let amin_str = args
        .next()
        .ok_or_else(|| KappaError::MissingArgument("--amin".to_string()))?;

    if amin_str.parse::<f64>().is_err() && Path::new(&amin_str).exists() {
        return Ok(SizeArg::File(amin_str));
    }

    let amin = amin_str
        .parse::<f64>()
        .map_err(|_| KappaError::InvalidType)?;

    let amax_str = args.next();
    if let Some(amax_str) = amax_str {
        if let Ok(amax) = amax_str.parse::<f64>() {
            let next_arg = args.next();
            if let Some(next_arg) = next_arg {
                if next_arg.contains(':') {
                    let parts: Vec<&str> = next_arg.split(':').collect();
                    if parts.len() != 2 {
                        return Err(KappaError::InvalidArgument(
                            "
                        Could not parse centroid size and logarithmic width"
                                .to_string(),
                        ));
                    }
                    let amean = parts[0]
                        .parse::<f64>()
                        .map_err(|_| KappaError::InvalidType)?;
                    let asigma = parts[1]
                        .parse::<f64>()
                        .map_err(|_| KappaError::InvalidType)?;
                    let na = args
                        .next()
                        .map(|s| s.parse::<usize>().map_err(|_| KappaError::InvalidType))
                        .transpose()?;
                    return Ok(SizeArg::LogNormal {
                        amin,
                        amax,
                        amean,
                        asigma,
                        na,
                    });
                } else if let Ok(apow) = next_arg.parse::<f64>() {
                    let na = args
                        .next()
                        .map(|s| s.parse::<usize>().map_err(|_| KappaError::InvalidType))
                        .transpose()?;
                    return Ok(SizeArg::PowerLaw {
                        amin,
                        amax: Some(amax),
                        apow: Some(apow),
                        na,
                    });
                } else {
                    return Err(KappaError::InvalidArgument(
                        "Could not parse powerlaw or lognormal arguments".to_string(),
                    ));
                }
            }
            return Ok(SizeArg::PowerLaw {
                amin,
                amax: Some(amax),
                apow: None,
                na: None,
            });
        } else {
            return Err(KappaError::InvalidType);
        }
    }

    Ok(SizeArg::PowerLaw {
        amin,
        amax: None,
        apow: None,
        na: None,
    })
}

fn process_wavelength() {
    todo!()
}

fn process_material<I>(
    material: &str,
    args: &mut Peekable<I>,
    materials: &mut Vec<String>,
    materials_mfrac: &mut Vec<f64>,
    material_types: &mut Vec<MaterialType>,
    material_type: MaterialType,
) -> Result<Option<(f64, f64, f64)>, KappaError>
where
    I: Iterator<Item = String>,
{
    if MATERIAL_KEYS.contains(&material) {
        materials.push(material.to_string());
        material_types.push(material_type);
        // See if there is a next argument
        // If there is, it should be a mass fraction
        if let Some(next_arg) = args.peek() {
            // If it is a number, it is a mass fraction
            if let Ok(mfrac) = next_arg.parse::<f64>() {
                materials_mfrac.push(mfrac);
                args.next();
                Ok(None)
            } else {
                // Next argument is not a number, so it is a material type
                materials_mfrac.push(1.0); // Default is to use 100% of the mass
                Ok(None)
            }
        } else {
            // This material has a default mass fraction of 100%
            materials_mfrac.push(1.0);
            // Move on to the next material
            Ok(None)
        }
    } else if material.ends_with(".lnk") {
        // Positional argument is a file path
        // Placeholder for now
        let rho_in = Option::<f64>::None;
        let _ = read_lnk_file(material, rho_in);
        Ok(None)
    } else if material.contains(':') {
        // This is n:k:rho format
        let parts: Vec<&str> = material.split(':').collect();
        if parts.len() != 3 {
            eprintln!("Invalid material format");
            return Err(KappaError::InvalidArgument("-c/-m".to_string()));
        } else {
            let n = parts[0]
                .parse::<f64>()
                .map_err(|_| KappaError::InvalidType)?;
            let k = parts[1]
                .parse::<f64>()
                .map_err(|_| KappaError::InvalidType)?;
            let rho = parts[2]
                .parse::<f64>()
                .map_err(|_| KappaError::InvalidType)?;
            materials.push(material.to_string());
            material_types.push(MaterialType::Core);
            materials_mfrac.push(1.0);
            Ok(Some((n, k, rho)))
        }
    } else {
        eprintln!("{} is not a valid material key", material);
        Err(KappaError::InvalidArgument("-c/-m".to_string()))
    }
}

fn print_short_help() {
    // Description
    print!("{}", "kappa".bold().color(Color::Yellow));
    print!(" : ");
    println!("Dust opacities from the command line");
    println!("        Dominik, Min, Tazaki 2021, https://ascl.net/2104.010, version 1.9.14");
    println!();

    // Usage
    print!("{}", "Usage".bold().underline().color(Color::Green));
    println!(":");
    println!("  kappa [OPTIONS] [MATERIAL [MFRAC]]...");
    println!();

    // Arguments
    print!("{}", "Arguments".bold().underline().color(Color::Cyan));
    println!(":");
    print!("{:<25}", "  <MATERIAL> [MFRAC]...");
    println!("Specify materials. Example: c-gra 0.5 c-z 0.3");
    println!();

    // Options
    print!("{}", "Options".bold().underline().color(Color::Blue));
    println!(":");

    // General options (grain composition)
    print!("{}", "GRAIN COMPOSITION".color(Color::BrightMagenta));
    println!(":");
    println!("{}", "  -L, --list".bold());
    println!("      List built-in materials");
    println!("{}", "  -m".bold().to_string() + " <MATERIAL> [MFRAC]...");
    println!(
        "      Specify a material to include in the mantle, with its mass fraction (default: 0.0)"
    );
    println!(
        "{}",
        "  -p".bold().to_string() + " <CORE_POROSITY> [MANTLE_POROSITY]..."
    );
    println!("      Specify core and mantle porosities (default: 0.0)");
    println!("{}", "  --diana".bold());
    println!("      Use the DIANA composition (Woitke et al. 2016)");
    println!("{}", "  --dsharp".bold());
    println!("      Use the DSHARP composition (Birnstiel et al. 2018)");
    println!("{}", "  --dsharp-no-ice".bold());
    println!("      Use the DSHARP composition without ice");
    println!();

    // Grain geometry and computational method options
    print!(
        "{}",
        "GRAIN GEOMETRY AND COMPUTATIONAL METHOD".color(Color::BrightMagenta)
    );
    println!(":");
    println!("{}", "  --dhs".bold().to_string() + " [FMAX]");
    println!("      Use the Distribution of Hollow Spheres (DHS) approach to model deviations");
    println!(
        "{}",
        "  --mmf".bold().to_string() + " [A0 [DFRAC-OR-FILL [KF]]]"
    );
    println!("      Use the Modified Mean Field (MMF) theory from Tazaki & Tanaka (2017)");
    println!("{}", "  --mie".bold());
    println!("      Do a standard Mie calculation for perfect spheres");
    println!("      Equivalent to --dhs 0.0");
    println!("{}", "  --cde".bold());
    println!(
        "      Compute Continuous Distribution of Ellipsoids (CDE) opacities in the Rayleigh limit"
    );
    println!();

    // Grain size distribution options
    print!("{}", "GRAIN SIZE DISTRIBUTION".color(Color::BrightMagenta));
    println!(":");
    println!(
        "{}",
        "  -a".bold().to_string() + " <AMIN> [AMAX [APOW [NA]]]"
    );
    println!("      Specify a power-law grain size distribution");
    println!(
        "{}",
        "  -a".bold().to_string() + " <AMIN AMAX AMEAN:ASIG> [NA]"
    );
    println!("      Specify a [log-]normal grain size distribution");
    println!("{}", "  -a".bold().to_string() + " <FILE>");
    println!("      Read grain size distribution from a file");
    println!("{}", "  --amin".bold().to_string() + " <AMIN>");
    println!("      Override minimum grain size");
    println!("{}", "  --amax".bold().to_string() + " <AMAX>");
    println!("      Override maximum grain size");
    println!("{}", "  --amean".bold().to_string() + " <AMEAN>");
    println!("      Override centroid size");
    println!("{}", "  --asig".bold().to_string() + " <ASIG>");
    println!("      Override logarithmic width");
    println!("{}", "  --apow".bold().to_string() + " <APOW>");
    println!("      Override power-law index");
    println!("{}", "  --na".bold().to_string() + " <NA>");
    println!("      Override number of size bins");
    println!();

    // Wavelength options
    print!("{}", "WAVELENGTH GRID".color(Color::BrightMagenta));
    println!(":");
    println!("{}", "  -l".bold().to_string() + " <LMIN> [LMAX [NLAM]]");
    println!("      Specify a wavelength grid");
    println!("{}", "  -l".bold().to_string() + " <FILE>");
    println!("      Read wavelength grid from a file");
    println!("{}", "  --lmin".bold().to_string() + " <LMIN>");
    println!("      Override minimum wavelength");
    println!("{}", "  --lmax".bold().to_string() + " <LMAX>");
    println!("      Override maximum wavelength");
    println!("{}", "  --nlam".bold().to_string() + " <NLAM>");
    println!("      Override number of wavelength points");
    println!();

    // Output options
    print!("{}", "OUTPUT".color(Color::BrightMagenta));
    println!(":");
    println!("{}", "  -o".bold().to_string() + " [DIR]");
    println!("      Write output to a directory (default: current directory)");
    println!("{}", "  -s".bold().to_string() + " [NANG]");
    println!("      Include scattering matrix in the output");
    println!("{}", "  -d".bold().to_string() + " [NSUB]");
    println!("      Divide computation upto into n_a parts");
    println!("{}", "  --chop".bold().to_string() + " [NDEG]");
    println!("      Cap the first NDEG degrees of the forward scattering peak");
    println!("{}", "  --fits".bold());
    println!("      Write output in FITS format");
    println!("{}", "  --radmc".bold().to_string() + " [LABEL]");
    println!("      Write output in RADMC-3D format");
    println!("{}", "  -w".bold());
    println!("      Create files on which kappa can operate");
    println!();

    // Miscellaneous options
    print!("{}", "MISCELLANEOUS".color(Color::BrightMagenta));
    println!(":");
    println!("{}", "  -h, --help".bold());
    println!("      Print this help message");
    println!("{}", "  -q, --quiet".bold());
    println!("      Quiet mode");
    println!("{}", "  -v, --version".bold());
    println!("      kappa version");
}

fn print_long_help() {
    todo!()
}

pub fn list_materials() {
    println!(
        ":===================== *** LIST OF BUILT-IN MATERIALS *** =====================:\n\
         : amorph.pyroxenes  pyr-mg100/95/80/70/60/50/40                                :\n\
         : amorph.olivines   ol-mg100/40                        (Dorschner95,Henning96) :\n\
         : cryst. pyr/ol     pyr-c-mg96 ol-c-mg100/95/00      (Jäger96,Suto06,Fabian01) :\n\
         : other silicates   astrosil                                        (Draine03) :\n\
         : amorphous carbon  c-z    c-p                           (Zubko96,Preibisch93) :\n\
         : graphite,special  c-gra  c-nano  c-org                (Dra.03,Mut.04,Hen.96) :\n\
         : quartz,corundum   sio2   cor-c                          (Kitamura07,Koike95) :\n\
         : iron/sulfide      fe-c   fes                                     (Henning96) :\n\
         : carbides          sic                                             (Draine93) :\n\
         : water ice         h2o-w  h2o-a                          (Warren08,Hudgins93) :\n\
         : other ices        co2-w  nh3-m                       (Warren86,Martonchik83) :\n\
         :                   co-a   co2-a/c ch4-a/c ch3oh-a/c   (Palumbo06,Gerakines20) :\n\
         :                                                                              :\n\
         : *** ABBREVIATIONS AND GENERIC KEYS FOR QUICK ACCESS ***                      :\n\
         : pyr  -> pyr-mg70     c    -> c-z         iron -> fe-c      h2o  -> h2o-w     :\n\
         : ol   -> ol-mg50      gra  -> c-gra       qua  -> sio2      co   -> co-a      :\n\
         : ens  -> pyr-c-mg96   org  -> c-org       cor  -> cor-c     co2  -> co2-w     :\n\
         : for  -> ol-c-mg100                       tro  -> fes       nh3  -> nh3-m     :\n\
         : fay  -> ol-c-mg00                                                            :\n\
         :==============================================================================:",
    );
}
