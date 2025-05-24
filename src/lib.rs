#![doc(
    html_logo_url = "https://raw.githubusercontent.com/ProsiaLAB/prosialab.github.io/refs/heads/main/images/prosialab.jpeg"
)]
//! # Introduction
//! This tool produces complex dust particle opacities right from the command line. It is
//! derived from Michiel Min’s DHS OpacityTool and also implements Ryo Tazaki’s MMF
//! theory for highly porous aggregates.
//!
//! For a brief overview behind the details of the implementation from a theoretical
//! perspective, read the [book](https://prosialab.github.io/kappa/book/index.html).
//!
//! ## Capabilities
//! - stand-alone tool, fully command line driven, no input files need to be edited
//! - full scattering matrix output in several formats, including for RADMC-3D
//! - combining materials through mixing into a complex grain with porosity
//! - built-in: a curated collection of materials for applications in astronomy
//! - external refractive index data can be used just as easily
//! - computational methods: (i) DHS (Distribution of Hollow Spheres) for irregular
//!   grains and low-porosity aggregates. Standard Mie theory for perfect spheres
//!   is available as a limiting case. (ii) MMF (Modified Mean Field) theory for
//!   high-porosity/fractal aggregates. (iii) CDE approximation in the Rayleigh limit.
//! - Python interface module for plotting and post-processing
//!
//! ## Physical units used by `kappa`
//! Due to conventions in our field, the input and output of `kappa` uses the following units
//!
//! - grain sizes and wavelengths are in microns
//! - mass densities of materials are in g/cm^3
//! - opacities are in cm^2/g
//! - scattering matrices are in sr^-1 or cm^-1 g^-1 sr^-1
//!
//! # Examples
//! A simple grain made only of the default pyroxene, for the default grain size distribution
//! ($a^{-3.5}$ powerlaw from 0.05 to 3000μm), on the default wavelength grid (0.05μm to 1cm).
//!
//! # Installation
//! The easiest way to install `kappa` is to first get the source code from GitHub:
//! ```bash
//! git clone https://github.com/ProsiaLAB/kappa.git
//! ```
//! Then, you can build the executable with `cargo`:
//! ```bash
//! cd kappa
//! cargo build --release
//! ```
//! This will create a binary called `kappa` in the `target/release` directory.
//!
//! # Usage
//! The command line interface is documented in the help message:
//! ```bash
//! kappa --help
//! ```
//!
//! ## Grain composition
//! If no composition is specified, the DIANA composition is used by default.
//!
//! ### Core material
//! Specify a material to include in the grain. MATERIAL can be the key for a builtin
//! material, the path to an lnk file, or colon-separated numbers n:k:rho3. MFRAC is
//! the mass fraction (default 1.0) of the material. You can give up to 20 materials
//! to build up the grain. Mass fractions do not have to add up to one, they will be
//! renormalized. All materials will be mixed together using the Bruggeman rule, and
//! vacuum can be added through the porosity. -c in front of each material is optional.
//!
//! ### Mantle material
//! Like -c, but place this material into the grain mantle. Multiple mantle materials
//! will be mixed using the Bruggeman rule, and then that mix will be added to the
//! core using the Maxwell-Garnett rule. The -m is not optional, it must be present.
//!
//! ### Porosity
//! Porosity, the volume fraction of vacuum, a number smaller than 1. The default is
//! 0. A single value will apply to both core and mantle, but a second value will be
//! specific for the mantle (and may be 0).
//!
//! ## Grain geometry and computational method
//! If no method is explicitly specified, the default is -dhs 0.8, i.e. DHS with fmax=0.8.
//!
//! ### Distribution of Hollow Spheres
//! Use the Distribution of Hollow Spheres (DHS, Min+ 2005) approach to model deviations from
//! perfect spherical symmetry and low-porosity aggregates. Spheres with
//! inner holes with volume fractions between 0 and fmax (default 0.8) are averaged to
//! mimic irregularities. fmax=0 means to use solid spheres (Mie theory), i.e. perfectly
//! regular grains. For backward compatibility, -fmax can be used instead of -dhs
//!
//! ### Modified Mean Field theory
//! Use Modified Mean Field theory (MMF, Tazaki & Tanaka 2018) to compute opacities of
//! highly porous or fractal aggregates. -c, -m, and -p determine the composition of
//! monomers with radius A0 (default 0.1μm). Particles will be aggregates
//! with a compact size given by the -a switch, giving rise to N = a3/a30 monomers.
//! DFRAC-OR-FILL specifies either the fractal dimension (if >1) or the volume filling
//! factor (if <1). The default is 0.2. KF may be used to change the default prefactor.
//!
//! ### Mie theory
//! Do a standard Mie calculation for perfect spheres. This is short for -dhs 0
//!
//! ### Continuum Distribution of Ellipsoids
//! Compute CDE (continuous distribution of ellipsoids) Rayleigh limit opacities
//!
//! ## Grain size distribution
//! Grain size distributions can be specified using a powerlaw, a (log-)normal or from
//! a file.
//!
//! ### Powerlaw size distribution
//! Specify (minimum) grain radius, and optionally maximum grain radius, the size
//! distribution powerlaw and the number of size bins. You may also use options to
//! set individual values with -amin, -amax, -apow, -na. The defaults are 0.05 μm,
//! 3000 μm, 3.5, and 15 per size decade with a fixed minimum of 5, respectively.
//!
//! ### (Log-)normal size distribution
//! Specify the centroid size and the logarithmic width for a log-normal size distribution.
//! You may also use -amean and -asig options to set these values. If ASIG is
//! negative, create a normal distribution with that width (in μm) around AMEAN.
//!
//! ## Wavelength grid
//! Specify the (minimum) wavelength, and optionally the maximum wavelength and
//! the number of wavelengths points for the construction of the wavelength grid. The
//! default values are 0.05 μm, 10000 μm, and 300, respectively. You may also use the
//! options -lmin, -lmax, and -nlam (or -nl) to set individual values.
//! > If only one wavelength is specified with -l, then λmax=λmin and nλ=1 are implied.
//!
//! Grid can also be read from a file.
//!
//! ## Controlling the output
//! The standard output is the file dustkappa.dat, with the opacities and the asymmetry
//! parameter g. The following options control and extend the output.
//!
//! -   Put the output files in directory DIR instead of the current working directory.
//!     ./output will be used if -o is present but DIR is not specified.
//! -   Include the scattering matrix in the output. NANG may optionally change the
//!     number of equally-spaced angular grid points to cover the range of angles between
//!     0 and 180 degrees. The default for NANG is 180 and should normally be just fine.
//! -   Divide the computation up into na parts to produce a file for each grain size. Each
//!     size will be an average over a range of NSUB grains around the real size.
//! -   Cap the first NDEG (2 if unspecified) degrees of the forward scattering peak.
//! -   Write dustkappa.fits instead of ASCII output. With -d, write na files.
//! -   RADMC-3D uses a different angular grid and scattering matrix normalization. File
//!     names will contain LABEL if specified and have the extension .inp.
//! -   Write to STDOUT instead of files. The default is to write λ, κabs, κsca, κext, and
//!     g. Many other outputs are possible, run optool -print ? for a full list. For
//!     readability, a header line may be printed to STDERR, but STDOUT gets only numbers
//!     which can be used in pipes and for redirection. You can use this to extract a single
//!     value, for example the 850μm extinction opacity of grains between 1 and 3mm:
//!     optool -a 1000 3000 -l 850 -print kext
//! -   Write the files optool_sd.dat and optool_lam.dat with the grain size distribution
//!     and the wavelength grid, respectively. Also, write optool_mix.lnk with the result
//!     of mixing refractive index data. Exit without doing a computation.
//!
//!
//! # Material properties
//! `kappa` needs refractive index data to work.  For your convenience, a
//! useful list of materials is compiled into `kappa`. You can also find
//! and use other data.
//! To access one of the built-in materials, specify the corresponding key
//! string like `pyr-mg70`. In each material class we have selected a
//! useful default, accessible with an even simpler generic key (for
//! example, `pyr` is an alias for `pyr-mg70`). Most of the built-in
//! refractive index datasets have a reasonably wide wavelength coverage -
//! the few exceptions are highlighted by bold-face numbers.  If a
//! material is being used outside of the measured region, `kappa` will
//! still function, using extrapolated optical properties.
//!
//! Even the limited number of materials we have selected to include with
//! `kappa` can be daunting. To get started with some kind of standard
//! opacity, we recommend to work with pyroxene `pyr`, carbon
//! `c`, and, at low temperatures, water ice `h2o` (Woitke+
//! 2016). If you need to account for sulfur, you may want to include
//! troilite `tro` (Birnstiel+ 2016).
//!
//!
//! | **Key**                   | **Material**                       | **State**                 |     $\rho \mathrm{(g/cm^3)}$             | $\lambda_{min}$ ($\mu{}m$) | $\lambda_{max}$ ($\mu{}m$) | **Reference**   
//! |---------------------------|------------------------------------|---------------------------|------------------------------------------|----------------------------|----------------------------|-------------------
//! | `pyr-mg100`               | $\mathrm{MgSiO_3}$                 | Amorphous                 |      2.71                                |         0.2                |         500                | [Dorschner+95](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D)
//! | `pyr-mg95`                | $\mathrm{Mg_{0.95}Fe_{0.05}SiO_3}$ | Amorphous                 |      2.74                                |         0.2                |         500                | [Dorschner+95](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D)  
//! | `pyr-mg80`                | $\mathrm{Mg_{0.8}Fe_{0.2}SiO_3}$   | Amorphous                 |       2.9                                |         0.2                |         500                | [Dorschner+95](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D)  
//! | `pyr-mg70` (`pyr`)        | $\mathrm{Mg_{0.7}Fe_{0.3}SiO_3}$   | Amorphous                 |      3.01                                |         0.2                |         500                | [Dorschner+95](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D)  
//! | `pyr-mg60`                | $\mathrm{Mg_{0.6}Fe_{0.4}SiO_3}$   | Amorphous                 |       3.1                                |         0.2                |         500                | [Dorschner+95](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D)  
//! | `pyr-mg50`                | $\mathrm{Mg_{0.5}Fe_{0.5}SiO_3}$   | Amorphous                 |       3.2                                |         0.2                |         500                | [Dorschner+95](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D)  
//! | `pyr-mg40`                | $\mathrm{Mg_{0.4}Fe_{0.6}SiO_3}$   | Amorphous                 |       3.3                                |         0.2                |         500                | [Dorschner+95](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D)  
//! | `pyr-c-mg96` (`ens`)      | $\mathrm{Mg_{0.96}Fe_{0.04}SiO3}$  | Crystalline               |       2.8                                |       2.0                  |        99                  | [Jäger+98](https://ui.adsabs.harvard.edu/abs/1998A%26A...339..904J)     
//! | `ol-mg50` (`ol`)          | $\mathrm{MgFeSiO_4}$               | Amorphous                 |      3.71                                |         0.2                |         500                | [Dorschner+95](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D)  
//! | `ol-mg40`                 | $\mathrm{Mg_{0.8}Fe_{1.2}SiO_4}$   | Amorphous                 |      3.71                                |         0.2                |         500                | [Dorschner+95](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D)  
//! | `ol-c-mg100` (`for`)      | $\mathrm{Mg_{2}SiO_4}$             | Crystalline               |      3.27                                |       5.0                  |         200                | [Suto+06](https://ui.adsabs.harvard.edu/abs/2006MNRAS.370.1599S)      
//! | `ol-c-mg95`               | $\mathrm{Mg_{1.9}Fe_{0.1}SiO_4}$   | Crystalline               |      3.33                                |       2.0                  |        8190                | [Fabian+01](https://ui.adsabs.harvard.edu/abs/2001A%26A...378..228F)     
//! | `ol-c-mg00` (`fay`)       | $\mathrm{Fe_{2}SiO_4}$             | Crystalline               |      4.39                                |       3.0                  |         250                | [Fabian+01](https://ui.adsabs.harvard.edu/abs/2001A%26A...378..228F)     
//! | [`astrosil`](ASTROSIL)    | $\mathrm{MgFeSiO_4}$               | Mixed                     |       3.3                                |        $6\times 10^-5$     |         $1 \times 10^5$    | [Draine+03](https://ui.adsabs.harvard.edu/abs/2003ApJ...598.1017D)
//! | `c-z` (`c`)               | $\mathrm{C}$                       | Amorphous?                |       1.8                                |        0.05                |         $1 \times 10^4$    | [Zubko+96](https://ui.adsabs.harvard.edu/abs/1996MNRAS.282.1321Z)      
//! | `c-p`                     | $\mathrm{C}$                       | Amorphous                 |       1.8                                |        0.11                |         800                | [Preibisch+93](https://ui.adsabs.harvard.edu/abs/1993A%26A...279..577P)
//! | `c-gra` (`gra`)           | $\mathrm{C}$ graphite              | Crystalline               |     2.16?                                |       0.001                |        1000                | [Draine+03](https://ui.adsabs.harvard.edu/abs/2003ApJ...598.1026D)     
//! | `c-org` (`org`)           | $\mathrm{CHON}$ organics           | Amorphous                 |       1.4                                |         0.1                |         $1 \times 10^5$    | [Henning+96](https://ui.adsabs.harvard.edu/abs/1996A%26A...311..291H)    
//! | `fe-c` (`iron`)           | $\mathrm{Fe}$                      | Metal                     |      7.87                                |         0.1                |         $1 \times 10^5$    | [Henning+96](https://ui.adsabs.harvard.edu/abs/1996A%26A...311..291H)     
//! | [`c_nano`](C_NANO)        | $\mathrm{C}$ nano-diamond          | Crystalline               |       2.3                                |        0.02                |       110                  | [Mutschke+04](https://ui.adsabs.harvard.edu/abs/2004A%26A...423..983M)   
//! | `fes` (`tro`)             | $\mathrm{FeS}$                     | Metal                     |      4.83                                |         0.1                |         $1 \times 10^5$    | [Henning+96](https://ui.adsabs.harvard.edu/abs/1996A%26A...311..291H)     
//! | `sic`                     | $\mathrm{SiC}$                     | Crystalline               |      3.22                                |       0.001                |        1000                | [Laor93](https://ui.adsabs.harvard.edu/abs/1993ApJ...402..441L)         
//! | `sio2` (`qua`)            | $\mathrm{SiO_2}$                   | Amorphous                 |      2.65                                |      0.0006                |         500                | [Kitamura+07](https://ui.adsabs.harvard.edu/abs/2007ApOpt..46.8118K)    
//! | `cor-c` (`cor`)           | $\mathrm{Al_{2}O_3}$               | Crystalline               |       4.0                                |         0.5                |        40                  | [Koike+95](https://ui.adsabs.harvard.edu/abs/1995Icar..114..203K)       
//! | `h2o-w` (`h2o`)           | Water ice                          | Crystalline               |      0.92                                |        0.04                |         $2 \times 10^6$    | [Warren+08](https://ui.adsabs.harvard.edu/abs/2008JGRD..11314220W)      
//! | `h2o-a`                   | Water ice                          | Amorphous                 |      0.92                                |        0.04                |         $2 \times 10^6$    | [Hudgins+93](https://ui.adsabs.harvard.edu/abs/1993ApJS...86..713H)     
//! | `co2-w` (`co2`)           | $\mathrm{CO_2}$ ice                | Crystalline               |       1.6                                |        0.05                |         $2 \times 10^5$    | [Warren+86](https://ui.adsabs.harvard.edu/abs/1986ApOpt..25.2650W)      
//! | `nh3-m` (`nh3`)           | $\mathrm{NH_3}$ ice                | Crystalline               |      0.75                                |        0.14                |         $2 \times 10^0$    | [Martonchik+83](https://ui.adsabs.harvard.edu/abs/1984ApOpt..23..541M)  
//! | `co-a` (`co`)             | $\mathrm{CO}$ ice                  | Amorphous                 |      0.81                                |       3.8                  |       5.8                  | [Palumbo+06](https://ui.adsabs.harvard.edu/abs/2006PCCP....8..279P)     
//! | `co2-a` / `c`             | $\mathrm{CO_2}$ ice                | Amorphous / Crystalline   |       1.2                                |       2.5                  |        20                  | [Gerakines+20](https://ui.adsabs.harvard.edu/abs/2020ApJ...901...52G)
//! | `ch4-a` / `c`             | $\mathrm{CH_4}$ ice                | Amorphous / Crystalline   |      0.47                                |       2.0                  |        20                  | [Gerakines+20](https://ui.adsabs.harvard.edu/abs/2020ApJ...901...52G)
//! | `ch3oh-a` / `c`           | $\mathrm{CH_{3}OH}$ ice            | Amorphous / Crystalline   | 0.78/1.02                                |       2.0                  |        24                  | [Gerakines+20](https://ui.adsabs.harvard.edu/abs/2020ApJ...901...52G)

pub mod cli;
pub mod components;
pub mod config;
pub mod dhs;
pub mod fractal;
pub mod geofractal;
pub mod io;
pub mod mie;
pub mod opac;
pub mod types;
pub mod utils;

#[allow(unused_imports)]
use crate::components::*;
