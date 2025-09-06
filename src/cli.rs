//! Command line interface for the `kappa`.
use std::env;
use std::iter::Peekable;
use std::path::Path;
use std::process::exit;

use anyhow::anyhow;
use anyhow::Result;
use colored::{Color, Colorize};
use extensions::types::RVector;

use crate::components::MATERIAL_KEYS;
use crate::io::read_wavelength_grid;
use crate::io::{read_lnk_file, read_sizedis_file};
use crate::opac::RefractiveIndexKind;
use crate::opac::{KappaConfig, KappaError, KappaMethod, SpecialConfigs};
use crate::opac::{Material, MaterialKind};
use crate::opac::{SizeDistribution, WavelengthKind};
use crate::utils::regrid_lnk_data;

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

enum WavelengthArg {
    Value {
        lmin: f64,
        lmax: Option<f64>,
        nlam: Option<usize>,
    },
    File(String),
}

pub fn launch() -> Result<KappaConfig, KappaError> {
    let mut kpc = KappaConfig::default();

    let mut args = env::args().skip(1).peekable();
    if args.peek().is_none() {
        // If no arguments are given, use the DIANA composition
        kpc = KappaConfig::diana();
        return Ok(kpc);
    }

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-c" => {
                if let Some(material_arg) = args.next() {
                    match process_material(
                        &kpc,
                        &material_arg.to_string(),
                        &mut args,
                        MaterialKind::Core,
                    ) {
                        Ok(material) => {
                            kpc.nmat += 1;
                            kpc.ncore += 1;
                            kpc.materials.push(material);
                        }
                        Err(KappaError::ZeroMassFraction) => {
                            args.next();
                        }
                        Err(e) => return Err(e),
                    }
                } else {
                    return Err(anyhow!(
                        "Missing material key/file/refractive index path after -c/-m"
                    )
                    .into());
                }
            }
            "-m" => {
                // Same as "-c"
                if let Some(material_arg) = args.next() {
                    match process_material(
                        &kpc,
                        &material_arg.to_string(),
                        &mut args,
                        MaterialKind::Mantle,
                    ) {
                        Ok(material) => {
                            kpc.nmat += 1;
                            kpc.nmant += 1;
                            kpc.materials.push(material);
                        }
                        Err(KappaError::ZeroMassFraction) => {
                            args.next();
                        }
                        Err(e) => return Err(e),
                    }
                } else {
                    return Err(anyhow!(
                        "Missing material key/file/refractive index path after -c/-m"
                    )
                    .into());
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
                    kpc.sizedis = SizeDistribution::Apow;
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
                    kpc.sizedis = SizeDistribution::Normal;
                }
                SizeArg::File(file) => {
                    // Read file
                    (kpc.na, kpc.ameans_file) = read_sizedis_file(&file)?;
                    kpc.sizedis = SizeDistribution::File;
                }
            },
            "--amin" => {
                kpc.amin = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --amin"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --amin"))?;
            }
            "--amax" => {
                kpc.amax = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --amax"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --amax"))?;
            }
            "--amean" => {
                kpc.amean = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --amean"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --amean"))?;
                kpc.sizedis = SizeDistribution::Normal;
            }
            "--asig" | "--asigma" => {
                kpc.asigma = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --asig"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --asig"))?;
                kpc.sizedis = SizeDistribution::Normal;
            }
            "--apow" => {
                kpc.apow = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --apow"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --apow"))?;
            }
            "--na" => {
                kpc.na = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --na"))?
                    .parse::<usize>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --na"))?;
            }
            // Wavelength options
            "-l" => match process_wavelength(&mut args)? {
                WavelengthArg::Value { lmin, lmax, nlam } => {
                    kpc.lmin = lmin;
                    if let Some(lmax) = lmax {
                        kpc.lmax = lmax;
                    }
                    if let Some(nlam) = nlam {
                        kpc.nlam = nlam;
                    }
                    kpc.wavelength_kind = WavelengthKind::Other;
                }
                WavelengthArg::File(file) => {
                    // Read file
                    let (nlam, lam) = read_wavelength_grid(&file)?;
                    kpc.nlam = nlam;
                    let min_lam = lam.iter().reduce(|a, b| if a < b { a } else { b });
                    let max_lam = lam.iter().reduce(|a, b| if a > b { a } else { b });
                    match (min_lam, max_lam) {
                        (Some(min), Some(max)) => {
                            kpc.lmin = *min;
                            kpc.lmax = *max;
                        }
                        _ => {
                            return Err(anyhow!("Could not read wavelength grid from file").into());
                        }
                    }
                    kpc.lam = lam;
                    kpc.wavelength_kind = WavelengthKind::File;
                }
            },
            "--lmin" => {
                kpc.lmin = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --lmin"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --lmin"))?;
            }
            "--lmax" => {
                kpc.lmax = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --lmax"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --lmax"))?;
            }
            "--nlam" => {
                kpc.nlam = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --nlam"))?
                    .parse::<usize>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --nlam"))?;
            }
            // Grain geometry and computational method options
            "-p" | "--porosity" => {
                kpc.pcore = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --porosity"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --porosity"))?;
                kpc.pmantle = match args.next() {
                    Some(arg) => arg.parse::<f64>().map_err(|_| {
                        anyhow!("Invalid value type for argument for mantle porosity")
                    })?,
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
                    .ok_or_else(|| anyhow!("Missing argument: --xlim"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --xlim"))?;
            }
            "--xlim-dhs" => {
                kpc.xlim_dhs = args
                    .next()
                    .ok_or_else(|| anyhow!("Missing argument: --xlim-dhs"))?
                    .parse::<f64>()
                    .map_err(|_| anyhow!("Invalid value type for argument: --xlim-dhs"))?;
            }
            "--mmf" | "--mmf-ss" => {
                kpc.method = KappaMethod::MMF;
                if arg == "--mmf-ss" {
                    kpc.mmf_ss = true;
                }
                let next_arg = args.next();
                if let Some(next_arg) = next_arg {
                    if let Ok(mmf_a0) = next_arg.parse::<f64>() {
                        kpc.mmf_a0 = mmf_a0;
                        if let Some(next_arg) = args.next() {
                            if let Ok(mmf_struct) = next_arg.parse::<f64>() {
                                kpc.mmf_struct = mmf_struct;
                                if let Some(next_arg) = args.next() {
                                    kpc.mmf_kf = next_arg
                                        .parse::<f64>()
                                        .map_err(|_| anyhow!("Unable to parse `mmf_kf`"))?;
                                } else {
                                    kpc.mmf_kf = 0.0;
                                }
                            } else {
                                return Err(anyhow!("Unable to parse `mmf_struct`").into());
                            }
                        } else {
                            kpc.mmf_struct = 0.2; // Default is a filling factor of 20%
                        }
                    } else {
                        return Err(anyhow!("Unable to parse `mmf_a0`").into());
                    }
                } else {
                    kpc.mmf_a0 = 0.1;
                }
            }
            "--cde" => {
                kpc.method = KappaMethod::CDE;
            }
            // Other options
            "-o" => {
                if let Some(next_arg) = args.next() {
                    kpc.outdir = next_arg.to_string();
                } else {
                    kpc.outdir = "output".to_string();
                }
            }
            "-s" | "--scat" | "--scatter" => {
                kpc.write_scatter = true;
                if let Some(nang) = args.next() {
                    kpc.nang = nang.parse::<usize>().map_err(|_| {
                        anyhow!("Unable to parse `nang` argument for -s/--scat/--scatter")
                    })?;
                }
            }
            "--sp" | "--sparse" => {
                kpc.write_scatter = true;
                if let Some(l1) = args.next() {
                    if let Ok(l1) = l1.parse::<f64>() {
                        if let Some(l2) = args.next() {
                            if let Ok(l2) = l2.parse::<f64>() {
                                kpc.nsparse += 1;
                                if kpc.nsparse > 10 {
                                    return Err(anyhow!("Too many sparse files").into());
                                }
                                if l1 < l2 {
                                    kpc.scatlammin[kpc.nsparse - 1] = l1;
                                    kpc.scatlammax[kpc.nsparse - 1] = l2;
                                } else {
                                    kpc.scatlammin[kpc.nsparse - 1] = l2;
                                    kpc.scatlammax[kpc.nsparse - 1] = l1;
                                }
                            } else {
                                return Err(anyhow!(
                                    "Parsing error in second argument of --sparse"
                                )
                                .into());
                            }
                        } else {
                            // One specific wavelength
                            kpc.nsparse += 1;
                            kpc.scatlammin[kpc.nsparse - 1] = l1;
                        }
                    } else {
                        return Err(anyhow!("Parsing error in first argument of --sparse").into());
                    }
                } else {
                    // No arguments
                    return Err(anyhow!("Missing arguments to --sparse").into());
                }
            }
            "--chop" => {
                if let Some(chop_angle) = args.next() {
                    kpc.chop_angle = chop_angle
                        .parse::<f64>()
                        .map_err(|_| anyhow!("Unable to parse `chop_angle` argument for --chop"))?;
                } else {
                    kpc.chop_angle = 2.0;
                }
            }
            "--radmc" | "--radmc3d" => {
                todo!()
            }
            "--fits" => {
                kpc.write_fits = true;
            }
            "-d" => {
                kpc.split = true;
                if let Some(nsub) = args.next() {
                    kpc.nsubgrains = nsub
                        .parse::<usize>()
                        .map_err(|_| anyhow!("Unable to parse `nsub` argument for -d/--split"))?;
                } else {
                    kpc.nsubgrains = 1;
                }
            }
            "-w" => {
                kpc.write_grid = true;
            }
            "-b" | "--blend-only" => {
                kpc.blend_only = true;
            }
            // Miscellaneous options
            "-L" | "--list" => {
                list_materials();
                return Ok(kpc);
            }
            "-h" => {
                print_short_help();
                exit(0);
            }
            "--help" => {
                print_long_help();
                exit(0);
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
            _ => match process_material(&kpc, arg.as_str(), &mut args, MaterialKind::Core) {
                Ok(material) => kpc.materials.push(material),
                Err(KappaError::ZeroMassFraction) => {
                    args.next();
                }
                Err(e) => return Err(e),
            },
        }
    }

    Ok(kpc)
}

fn process_grain_size<I>(args: &mut I) -> Result<SizeArg, KappaError>
where
    I: Iterator<Item = String>,
{
    let amin_str = args
        .next()
        .ok_or_else(|| anyhow!("Missing argument: --amin"))?;

    if amin_str.parse::<f64>().is_err() && Path::new(&amin_str).exists() {
        return Ok(SizeArg::File(amin_str));
    }

    let mut amin = amin_str
        .parse::<f64>()
        .map_err(|_| anyhow!("Invalid value type for argument: --amin"))?;

    let amax_str = args.next();
    if let Some(amax_str) = amax_str {
        if let Ok(mut amax) = amax_str.parse::<f64>() {
            if amax < 0.0 {
                if amin + amax <= 0.0 {
                    return Err(anyhow!("Invalid size range: amin + amax <= 0").into());
                }
                amin += amax;
                amax = amin - 2.0 * amax;
                let apow = 0.0;
                let na: Option<usize> = None;
                return Ok(SizeArg::PowerLaw {
                    amin,
                    amax: Some(amax),
                    apow: Some(apow),
                    na,
                });
            }
            let next_arg = args.next();
            if let Some(next_arg) = next_arg {
                if next_arg.contains(':') {
                    let parts: Vec<&str> = next_arg.split(':').collect();
                    if parts.len() != 2 {
                        return Err(
                            anyhow!("Could not parse centroid size and logarithmic width").into(),
                        );
                    }
                    let amean = parts[0]
                        .parse::<f64>()
                        .map_err(|_| anyhow!("Invalid value type for argument `amean`"))?;
                    let asigma = parts[1]
                        .parse::<f64>()
                        .map_err(|_| anyhow!("Invalid value type for argument `asigma`"))?;
                    let na = args
                        .next()
                        .map(|s| {
                            s.parse::<usize>()
                                .map_err(|_| anyhow!("Invalid value type for argument `na`"))
                        })
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
                        .map(|s| {
                            s.parse::<usize>()
                                .map_err(|_| anyhow!("Invalid value type for argument `na`"))
                        })
                        .transpose()?;
                    return Ok(SizeArg::PowerLaw {
                        amin,
                        amax: Some(amax),
                        apow: Some(apow),
                        na,
                    });
                } else {
                    return Err(anyhow!("Could not parse powerlaw or lognormal arguments").into());
                }
            }
            Ok(SizeArg::PowerLaw {
                amin,
                amax: Some(amax),
                apow: None,
                na: None,
            })
        } else {
            Err(anyhow!("Error in parsing amax").into())
        }
    } else {
        let amax = amin;
        let na = 1;
        Ok(SizeArg::PowerLaw {
            amin,
            amax: Some(amax),
            apow: None,
            na: Some(na),
        })
    }
}

fn process_wavelength<I>(args: &mut I) -> Result<WavelengthArg, KappaError>
where
    I: Iterator<Item = String>,
{
    let lmin_str = args
        .next()
        .ok_or_else(|| anyhow!("Missing argument: --lmin"))?;

    if lmin_str.parse::<f64>().is_err() && Path::new(&lmin_str).exists() {
        return Ok(WavelengthArg::File(lmin_str));
    }

    let lmin = lmin_str
        .parse::<f64>()
        .map_err(|_| anyhow!("Invalid value type for argument: --lmin"))?;

    let lmax_str = args.next();
    if let Some(lmax_str) = lmax_str {
        if let Ok(lmax) = lmax_str.parse::<f64>() {
            let next_arg = args.next();
            if let Some(next_arg) = next_arg {
                if let Ok(nlam) = next_arg.parse::<usize>() {
                    Ok(WavelengthArg::Value {
                        lmin,
                        lmax: Some(lmax),
                        nlam: Some(nlam),
                    })
                } else {
                    Err(anyhow!("Could not parse nlam").into())
                }
            } else {
                let nlam = 300;
                Ok(WavelengthArg::Value {
                    lmin,
                    lmax: Some(lmax),
                    nlam: Some(nlam),
                })
            }
        } else {
            Err(anyhow!("Could not parse lmax").into())
        }
    } else {
        let lmax = lmin;
        let nlam = 1;
        Ok(WavelengthArg::Value {
            lmin,
            lmax: Some(lmax),
            nlam: Some(nlam),
        })
    }
}

fn process_material<I>(
    kpc: &KappaConfig,
    material_arg: &str,
    args: &mut Peekable<I>,
    material_type: MaterialKind,
) -> Result<Material, KappaError>
where
    I: Iterator<Item = String>,
{
    let mut material = Material::default();
    if MATERIAL_KEYS.contains(&material_arg) {
        material.key = material_arg.to_string();
        material.kind = material_type;
        // let component = get_lnk_data(&material.key);
        // (material.re, material.im) =
        //     regrid_lnk_data(component.l0, component.n0, component.k0, &kpc.lam, true);
        // See if there is a next argument
        // If there is, it should be a mass fraction
        if let Some(next_arg) = args.peek() {
            // If it is a number, it is a mass fraction
            if let Ok(mfrac) = next_arg.parse::<f64>() {
                material.mfrac = if mfrac == 0.0 {
                    return Err(KappaError::ZeroMassFraction);
                } else {
                    mfrac
                };
                args.next();
                Ok(material)
            } else {
                // Next argument is not a number, so it is a material type
                material.mfrac = 1.0; // Default is to use 100% of the mass
                Ok(material)
            }
        } else {
            // This material has a default mass fraction of 100%
            material.mfrac = 1.0;
            // Move on to the next material
            Ok(material)
        }
    } else if material_arg.ends_with(".lnk") {
        // Positional argument is a file path
        // Placeholder for now
        let rho_in = Option::<f64>::None;
        let component = read_lnk_file(material_arg, rho_in)?;
        let l0_slice: &[f64] = &component.l0.to_vec();
        let n0_slice: &[f64] = &component.n0.to_vec();
        let k0_slice: &[f64] = &component.k0.to_vec();
        (material.re, material.im) = regrid_lnk_data(l0_slice, n0_slice, k0_slice, &kpc.lam, true);
        material.rho = component.rho;
        material.cmd = RefractiveIndexKind::File;
        Ok(material)
    } else if material_arg.contains(':') {
        // This is n:k:rho format
        material.cmd = RefractiveIndexKind::CmdLine;
        let parts: Vec<&str> = material_arg.split(':').collect();
        if parts.len() != 3 {
            Err(anyhow!("Invalid material format").into())
        } else {
            material.n = parts[0]
                .parse::<f64>()
                .map_err(|_| anyhow!("Invalid value type for argument `n`"))?;
            material.k = parts[1]
                .parse::<f64>()
                .map_err(|_| anyhow!("Invalid value type for argument `k`"))?;
            material.rho = parts[2]
                .parse::<f64>()
                .map_err(|_| anyhow!("Invalid value type for argument `rho`"))?;
            if material.rho <= 0.0 {
                eprintln!("Density must be positive");
                return Err(anyhow!("Density must be positive").into());
            }
            material.kind = material_type;
            material.re = material.n * RVector::ones(kpc.nlam);
            material.im = material.k * RVector::ones(kpc.nlam);
            Ok(material)
        }
    } else {
        Err(anyhow!("Invalid material key").into())
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

    println!("{}", "  -c".bold().to_string() + " <MATERIAL> [MFRAC]...");
    println!("      Specify a material to include in the core, with its mass fraction");
    println!("      (default: 1.0)");

    println!("{}", "  -m".bold().to_string() + " <MATERIAL> [MFRAC]...");
    println!("      Specify a material to include in the mantle, with its mass fraction");
    println!("      (default: 1.0)");
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
    println!("{}", "  -L, --list".bold());
    println!("      List built-in materials");
    println!("{}", "  -h, --help".bold());
    println!("      Print this help message");
    println!("{}", "  -q, --quiet".bold());
    println!("      Quiet mode");
    println!("{}", "  -v, --version".bold());
    println!("      kappa version");
}

fn print_long_help() {
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
    println!("{:<25}", "  <MATERIAL> [MFRAC]...");
    println!();
    println!("Specify a material to include in the grain. MATERIAL can be the key for a builtin");
    println!("material, the path to an lnk file, or colon-separated numbers `n:k:rho`. MFRAC is");
    println!("the mass fraction (default 1.0) of the material. You can give up to 20 materials");
    println!("to build up the grain. Mass fractions do not have to add up to one, they will be");
    println!("renormalized. All materials will be mixed together using the Bruggeman rule, and");
    println!("vacuum can be added through the porosity.");
    println!();
    println!("`-c` in front of each material is optional.");

    println!();

    // Options
    print!("{}", "Options".bold().underline().color(Color::Blue));
    println!(":");

    // General options (grain composition)
    print!("{}", "GRAIN COMPOSITION".color(Color::BrightMagenta));
    println!(":");

    println!("{}", "  -c".bold().to_string() + " <MATERIAL> [MFRAC]...");
    println!("      Specify a material to include in the core, with its mass fraction");
    println!();
    println!("{}", "      (default: 1.0)".italic());
    println!();
    println!("{}", "  -m".bold().to_string() + " <MATERIAL> [MFRAC]...");
    println!("      Specify a material to include in the mantle, with its mass fraction");
    println!("{}", "      (default: 1.0)".italic());
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
    println!("{}", "  -L, --list".bold());
    println!("      List built-in materials");
    println!("{}", "  -h, --help".bold());
    println!("      Print this help message");
    println!("{}", "  -q, --quiet".bold());
    println!("      Quiet mode");
    println!("{}", "  -v, --version".bold());
    println!("      kappa version");
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
