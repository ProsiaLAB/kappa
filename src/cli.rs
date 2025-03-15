//! Command line interface for the `kappa`.
use std::{collections::HashMap, env, iter::Peekable};

use anyhow::Result;
use colored::{Color, Colorize};

// use crate::io::read_lnk_file;
use crate::components::StaticComponent;
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

pub fn launch() -> Result<KappaConfig, KappaError> {
    let mut kpc = KappaConfig::default();

    let mut args = env::args().skip(1).peekable();
    if args.peek().is_none() {
        return Ok(());
    }

    let mut materials: HashMap<String, StaticComponent> = HashMap::new();

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-m" => {
                todo!()
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
            "-a" => {
                todo!()
            }
            "--amin" => {
                kpc.amin = args
                    .next()
                    .unwrap()
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidArgument)?;
            }
            "--amax" => {
                kpc.amax = args
                    .next()
                    .unwrap()
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidArgument)?;
            }
            "--amean" => {
                kpc.amean = args
                    .next()
                    .unwrap()
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidArgument)?;
            }
            "--asig" => {
                kpc.asigma = args
                    .next()
                    .unwrap()
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidArgument)?;
            }
            "--apow" => {
                kpc.apow = args
                    .next()
                    .unwrap()
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidArgument)?;
            }
            "--na" => {
                kpc.na = args
                    .next()
                    .unwrap()
                    .parse::<usize>()
                    .map_err(|_| KappaError::InvalidArgument)?;
            }
            // Wavelength options
            "-l" => {
                todo!()
            }
            "--lmin" => {
                kpc.lmin = args
                    .next()
                    .unwrap()
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidArgument)?;
            }
            "--lmax" => {
                kpc.lmax = args
                    .next()
                    .unwrap()
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidArgument)?;
            }
            "--nlam" => {
                kpc.nlam = args
                    .next()
                    .unwrap()
                    .parse::<usize>()
                    .map_err(|_| KappaError::InvalidArgument)?;
            }
            // Grain geometry and computational method options
            "-p" => {
                todo!()
            }
            "--dhs" | "--fmax" | "--mie" => {
                kpc.method = KappaMethod::DHS;
                if arg == "--mie" {
                    kpc.fmax = 0.0;
                } else if args.peek().is_some() {
                    kpc.fmax = args
                        .next()
                        .unwrap()
                        .parse::<f64>()
                        .map_err(|_| KappaError::InvalidArgument)?;
                } else {
                    kpc.fmax = 0.8;
                }
            }
            "--xlim" => {
                let xlim = args
                    .next()
                    .ok_or_else(|| KappaError::InvalidArgument("--xlim".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--xlim-dhs" => {
                let xlim_dhs = args
                    .next()
                    .ok_or_else(|| KappaError::InvalidArgument("--xlim-dhs".to_string()))?
                    .parse::<f64>()
                    .map_err(|_| KappaError::InvalidType)?;
            }
            "--mmf" | "--mmf-ss" => {
                kpc.method = KappaMethod::MMF;
                if arg == "--mmf-ss" {
                    let mmf_ss = args
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
                return Ok(());
            }
            "-h" => {
                print_short_help();
                return Ok(());
            }
            "--help" => {
                print_long_help();
                return Ok(());
            }
            "-q" | "--quiet" => {
                todo!()
            }
            "-v" | "--version" => {
                todo!()
            }
            // Unknown argument
            _ => {
                println!("Unknown argument: {}", arg);
                return Ok(());
            }
        }
    }
    // println!("args: {:?}", args);
    // todo!()
    Ok(())
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
    println!("Specify materials. Example: c-gra 0.5 c-sil 0.3");
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
