//! Command line interface for the `kappa`.
use std::env;

use anyhow::Result;
use colored::{Color, Colorize};

// use crate::io::read_lnk_file;
use crate::opac::KappaError;

pub fn launch() -> Result<(), KappaError> {
    let mut args = env::args().skip(1).peekable();
    if args.peek().is_none() {
        print_usage();
        return Ok(());
    }

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print_usage();
                return Ok(());
            }
            "-q" => {
                todo!()
            }
            "-v" => {
                todo!()
            }
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

fn print_usage() {
    print!("{}", "kappa".bold().color(Color::Yellow));
    print!(" : ");
    println!("Dust opacities from the command line");
    println!("Dominik, Min, Tazaki 2021, https://ascl.net/2104.010, version 1.9.14");
    println!();
    print!("{}", "Usage".bold().underline().color(Color::Green));
    println!(":");
    println!("kappa [OPTIONS] [MATERIAL [MFRAC]]...");
    println!();
    print!("{}", "Arguments".bold().underline().color(Color::Cyan));
    println!(":");
    print!("{:<25}", "MATERIAL [MFRAC]...");
    println!("Specify materials. Example: c-gra 0.5 c-sil 0.3");
    println!();
    print!("{}", "Options".bold().underline().color(Color::Blue));
    println!(":");
    print!("{:<40}", "-L, --list".bold());
    println!("List built-in materials");
    print!("{:<2}{:<33}", "-c".bold(), " MATERIAL [MFRAC]...");
    println!("Specify a material to include in the grain, with its mass fraction (default: 1.0)");
    print!("{:<2}{:<33}", "-m".bold(), " MATERIAL [MFRAC]...");
    println!("Specify a material to include in the mantle, with its mass fraction (default: 0.0)");
    print!(
        "{:<2}{:<33}",
        "-p".bold(),
        " CORE_POROSITY [MANTLE_POROSITY]"
    );
    println!("Specify core and mantle porosities (default: 0.0)");
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
