//! Command line interface for the `kappa`.

use anyhow::Result;
use clap::ArgAction;
use clap::{Arg, Command};

use crate::io::read_lnk_file;
use crate::opac::KappaError;

pub fn run() -> Result<(), KappaError> {
    let matches = Command::new("kappa")
        .version("0.1.0")
        .about(
            "Dust opacities from the command line\n\
                Dominik, Min, Tazaki 2021, https://ascl.net/2104.010, version 1.9.14",
        )
        .arg(
            Arg::new("list-materials")
                .short('L')
                .long("list")
                .num_args(0)
                .help("List available materials"),
        )
        .arg(
            Arg::new("add-material-one-by-one")
                .short('c')
                .num_args(1..=2) // Accepts one or more arguments
                .action(ArgAction::Append) // Allows multiple uses of `-c`
                .value_names(["MATERIAL", "[MFRAC]"])
                .help("Specify a material and its mass fractions (optionally)."),
        )
        .arg(
            Arg::new("materials")
                .num_args(1..)
                .action(ArgAction::Append)
                .help("Specify materials. Example: c-gra 0.5 c-sil 0.3"),
        )
        .arg(
            Arg::new("add-mantle-material")
                .short('m')
                .num_args(1..=2)
                .value_names(["MATERIAL", "[MFRAC]"])
                .help("Add material with mass fraction in mantle"),
        )
        .arg(
            Arg::new("porosity")
                .short('p')
                .num_args(1..)
                .help("Set porosity, possibly different for core and mantle"),
        )
        .arg(
            Arg::new("dhs")
                .short('H')
                .long("dhs")
                .num_args(1)
                .help("Maximum volume fraction of vacuum in DHS computation"),
        )
        .arg(
            Arg::new("mmf")
                .short('M')
                .long("mmf")
                .num_args(1..)
                .help("Use MMF with monomer size A0 and fractal dim or fill"),
        )
        .arg(
            Arg::new("grain-size")
                .short('a')
                .num_args(1..) // Accepts at least one argument (file or numbers)
                .help("Set grain size distribution or read from a file"),
        )
        .arg(
            Arg::new("wavelength-grid")
                .short('l')
                .num_args(1..)
                .help("Set up wavelength grid (units for a and l: microns)"),
        )
        .arg(
            Arg::new("write-na")
                .short('d')
                .num_args(1)
                .help("Write NA files for specific grain sizes"),
        )
        .arg(
            Arg::new("scattering-matrix")
                .short('s')
                .num_args(1)
                .help("Add scattering matrix for NANG angles to output"),
        )
        .arg(
            Arg::new("chop")
                .short('C')
                .long("chop")
                .num_args(1)
                .help("Remove NDEG degrees from the forward-scattering peak"),
        )
        .arg(
            Arg::new("output-directory")
                .short('o')
                .num_args(1)
                .help("Output to DIRECTORY instead of current working dir"),
        )
        .arg(
            Arg::new("radmc")
                .short('R')
                .long("radmc")
                .num_args(1)
                .help("Output as RADMC-3D input file"),
        )
        .arg(
            Arg::new("fits")
                .short('F')
                .long("fits")
                .action(ArgAction::SetTrue)
                .help("Output to FITS file"),
        )
        .arg(
            Arg::new("print")
                .short('P')
                .long("print")
                .num_args(1)
                .help("Output to STDOUT"),
        )
        .arg(
            Arg::new("quiet")
                .short('q')
                .long("quiet")
                .action(ArgAction::SetTrue)
                .help("Quiet mode"),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue)
                .help("More verbose on STDOUT"),
        )
        .arg(
            Arg::new("xlim")
                .short('X')
                .long("xlim")
                .num_args(1)
                .help("Switch to Mie for speed at x=XLIM"),
        )
        .get_matches();

    // check if first argument is a toml file
    // println!("Arguments: {:?}", matches.);
    // exit(1);

    // if let Some(materials) = matches.get_many::<String>("materials") {
    //     let materials: Vec<&String> = materials.collect();
    //     println!("Materials: {:?}", materials);
    // }

    // if let Some(materials) = matches.get_many::<String>("material") {
    //     let materials: Vec<&String> = materials.collect();
    //     println!("Materials: {:?}", materials);
    // }

    // Check if `-a` was provided with values
    if let Some(values) = matches.get_many::<String>("grain-size") {
        let values: Vec<String> = values.map(|s| s.to_string()).collect();

        if values.len() == 1 {
            // Make sure the value is a file
            if values[0].contains(".lnk") {
                println!("Reading size distribution from file: {}", values[0]);
                // Read the file
                _ = read_lnk_file(&values[0]).map_err(|_| KappaError::InvalidLNKFile)?;
            } else {
                return Err(KappaError::InvalidSizeInput);
            }
        } else if values.len() >= 2 && values[2].contains(':') {
            println!(
                "Using AMIN={}, AMAX={}, AMEAN:ASIG={} (Gaussian-like distribution)",
                values[0], values[1], values[2]
            );
            if values.len() > 3 {
                println!("Additional NA parameter: {}", values[3]);
            }
        } else {
            println!(
                "Using AMIN={}, AMAX={}, APOW={}, NA={}",
                values[0],
                values.get(1).map(|s| s.as_str()).unwrap_or("-"),
                values.get(2).map(|s| s.as_str()).unwrap_or("-"),
                values.get(3).map(|s| s.as_str()).unwrap_or("-"),
            );
        }
    }

    if matches.get_flag("list-materials") {
        list_materials();
    }

    Ok(())
}

fn list_materials() {
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
