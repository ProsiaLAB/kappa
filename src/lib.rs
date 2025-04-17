//! This tool produces complex dust particle opacities right from the command line. It is
//! derived from Michiel Min’s DHS OpacityTool and also implements Ryo Tazaki’s MMF
//! theory for highly porous aggregates.
//!
//! # Capabilities
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
//! # Physical units used by `kappa`
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
//! \fbox{c}, and, at low temperatures, water ice `h2o` (Woitke+
//! 2016). If you need to account for sulfur, you may want to include
//! troilite `tro` (Birnstiel+ 2016).
//!
//!
//! | *-c Key*   | *-c Key*    | *Material*              | *State*     |      \rho | \lambda_min | \lambda_max | *Reference*   | *Comment*    | *File*                      |
//! | generic    | full key    |                         |             |    g/cm^3 |      \mu{}m |      \mu{}m |               |              |                             |
//! |------------|-------------|-------------------------|-------------|-----------|-------------|-------------|---------------|--------------|-----------------------------|
//! |            | pyr-mg100   | MgSiO_3                 | amorph      |      2.71 |         0.2 |         500 | [[https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D][Dorschner+95]]  |              | [[file:lnk_data/pyr-mg100-Dorschner1995.lnk][pyr-mg100-Dorschner1995.lnk]] |
//! |            | pyr-mg95    | Mg_{0.95}Fe_{0.05}SiO_3 | amorph      |      2.74 |         0.2 |         500 | [[https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D][Dorschner+95]]  |              | [[file:lnk_data/pyr-mg95-Dorschner1995.lnk][pyr-mg95-Dorschner1995.lnk]]  |
//! |            | pyr-mg80    | Mg_{0.8}Fe_{0.2}SiO_3   | amorph      |       2.9 |         0.2 |         500 | [[https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D][Dorschner+95]]  | \rho interp. | [[file:lnk_data/pyr-mg80-Dorschner1995.lnk][pyr-mg80-Dorschner1995.lnk]]  |
//! | \fbox{pyr} | pyr-mg70    | Mg_{0.7}Fe_{0.3}SiO_3   | amorph      |      3.01 |         0.2 |         500 | [[https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D][Dorschner+95]]  |              | [[file:lnk_data/pyr-mg70-Dorschner1995.lnk][pyr-mg70-Dorschner1995.lnk]]  |
//! |            | pyr-mg60    | Mg_{0.6}Fe_{0.4}SiO_3   | amorph      |       3.1 |         0.2 |         500 | [[https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D][Dorschner+95]]  | \rho interp. | [[file:lnk_data/pyr-mg60-Dorschner1995.lnk][pyr-mg60-Dorschner1995.lnk]]  |
//! |            | pyr-mg50    | Mg_{0.5}Fe_{0.5}SiO_3   | amorph      |       3.2 |         0.2 |         500 | [[https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D][Dorschner+95]]  |              | [[file:lnk_data/pyr-mg50-Dorschner1995.lnk][pyr-mg50-Dorschner1995.lnk]]  |
//! |            | pyr-mg40    | Mg_{0.4}Fe_{0.6}SiO_3   | amorph      |       3.3 |         0.2 |         500 | [[https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D][Dorschner+95]]  | \rho interp. | [[file:lnk_data/pyr-mg40-Dorschner1995.lnk][pyr-mg40-Dorschner1995.lnk]]  |
//! | ens        | pyr-c-mg96  | Mg_{0.96}Fe_{0.04}SiO3  | cryst[fn:4] |       2.8 |       *2.0* |        *99* | [[https://ui.adsabs.harvard.edu/abs/1998A%26A...339..904J][Jäger+98]]      |              | [[file:lnk_data/pyr-c-mg96-Jäger1998.lnk][pyr-c-mg96-Jäger1998.lnk]]    |
//! |------------+-------------+-------------------------+-------------+-----------+-------------+-------------+---------------+--------------+-----------------------------|
//! | ol         | ol-mg50     | MgFeSiO_4               | amorph      |      3.71 |         0.2 |         500 | [[https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D][Dorschner+95]]  |              | [[file:lnk_data/ol-mg50-Dorschner1995.lnk][ol-mg50-Dorschner1995.lnk]]   |
//! |            | ol-mg40     | Mg_{0.8}Fe_{1.2}SiO_4   | amorph      |      3.71 |         0.2 |         500 | [[https://ui.adsabs.harvard.edu/abs/1995A%26A...300..503D][Dorschner+95]]  | \rho ?       | [[file:lnk_data/ol-mg40-Dorschner1995.lnk][ol-mg40-Dorschner1995.lnk]]   |
//! | for        | ol-c-mg100  | Mg_{2}SiO_4             | cryst[fn:4] |      3.27 |       *5.0* |         200 | [[https://ui.adsabs.harvard.edu/abs/2006MNRAS.370.1599S][Suto+06]]       | switch out?  | [[file:lnk_data/ol-c-mg100-Suto2006.lnk][ol-c-mg100-Suto2006.lnk]]     |
//! |            | ol-c-mg95   | Mg_{1.9}Fe_{0.1}SiO_4   | cryst[fn:4] |      3.33 |       *2.0* |        8190 | [[https://ui.adsabs.harvard.edu/abs/2001A%26A...378..228F][Fabian+01]]     | \rho ?       | [[file:lnk_data/ol-c-mg95-Fabian2001.lnk][ol-c-mg95-Fabian2001.lnk]]    |
//! | fay        | ol-c-mg00   | Fe_{2}SiO_4             | cryst[fn:4] |      4.39 |       *3.0* |         250 | [[https://ui.adsabs.harvard.edu/abs/2001A%26A...378..228F][Fabian+01]]     |              | [[file:lnk_data/ol-c-mg00-Fabian2001.lnk][ol-c-mg00-Fabian2001.lnk]]    |
//! |------------+-------------+-------------------------+-------------+-----------+-------------+-------------+---------------+--------------+-----------------------------|
//! |            | astrosil    | MgFeSiO_4               | mixed       |       3.3 |        6e-5 |         1e5 | [[https://ui.adsabs.harvard.edu/abs/2003ApJ...598.1017D][Draine+03]]     |              | [[file:lnk_data/astrosil-Draine2003.lnk][astrosil-Draine2003.lnk]]     |
//! |------------+-------------+-------------------------+-------------+-----------+-------------+-------------+---------------+--------------+-----------------------------|
//! | \fbox{c}   | c-z         | C                       | amorph?     |       1.8 |        0.05 |         1e4 | [[https://ui.adsabs.harvard.edu/abs/1996MNRAS.282.1321Z][Zubko+96]]      |              | [[file:lnk_data/c-z-Zubko1996.lnk][c-z-Zubko1996.lnk]]           |
//! |            | c-p         | C                       | amorph      |       1.8 |        0.11 |         800 | [[https://ui.adsabs.harvard.edu/abs/1993A%26A...279..577P][Preibisch+93]]  |              | [[file:lnk_data/c-p-Preibisch1993.lnk][c-p-Preibisch1993.lnk]]       |
//! | gra        | c-gra       | C graphite              | cryst[fn:4] |     2.16? |       0.001 |        1000 | [[https://ui.adsabs.harvard.edu/abs/2003ApJ...598.1026D][Draine+03]]     |              | [[file:lnk_data/c-gra-Draine2003.lnk][c-gra-Draine2003.lnk]]        |
//! | org        | c-org       | CHON organics           | amorph      |       1.4 |         0.1 |         1e5 | [[https://ui.adsabs.harvard.edu/abs/1996A%26A...311..291H][Henning+96]]    |              | [[file:lnk_data/c-org-Henning1996.lnk][c-org-Henning1996.lnk]]       |
//! |            | c-nano      | C nano-diamond          | cryst       |       2.3 |        0.02 |       *110* | [[https://ui.adsabs.harvard.edu/abs/2004A%26A...423..983M][Mutschke+04]]   |              | [[file:lnk_data/c-nano-Mutschke2004.lnk][c-nano-Mutschke2004.lnk]]     |
//! |------------+-------------+-------------------------+-------------+-----------+-------------+-------------+---------------+--------------+-----------------------------|
//! | iron       | fe-c        | Fe                      | metal       |      7.87 |         0.1 |         1e5 | [[https://ui.adsabs.harvard.edu/abs/1996A%26A...311..291H][Henning+96]]    |              | [[file:lnk_data/fe-c-Henning1996.lnk][fe-c-Henning1996.lnk]]        |
//! | \fbox{tro} | fes         | FeS                     | metal       |      4.83 |         0.1 |         1e5 | [[https://ui.adsabs.harvard.edu/abs/1996A%26A...311..291H][Henning+96]]    |              | [[file:lnk_data/fes-Henning1996.lnk][fes-Henning1996.lnk]]         |
//! |            | sic         | SiC                     | cryst       |      3.22 |       0.001 |        1000 | [[https://ui.adsabs.harvard.edu/abs/1993ApJ...402..441L][Laor93]]        |              | [[file:lnk_data/sic-Draine1993.lnk][sic-Draine1993.lnk]]          |
//! |------------+-------------+-------------------------+-------------+-----------+-------------+-------------+---------------+--------------+-----------------------------|
//! | qua        | sio2        | SiO_2                   | amorph      |      2.65 |      0.0006 |         500 | [[https://ui.adsabs.harvard.edu/abs/2007ApOpt..46.8118K][Kitamura+07]]   | \rho ?       | [[file:lnk_data/sio2-Kitamura2007.lnk][si02-Kitamura2007.lnk]]       |
//! | cor        | cor-c       | Al_{2}O_3               | cryst       |       4.0 |         0.5 |        *40* | [[https://ui.adsabs.harvard.edu/abs/1995Icar..114..203K][Koike+95]]      |              | [[file:lnk_data/cor-c-Koike1995.lnk][cor-c-Koike1995.lnk]]         |
//! |------------+-------------+-------------------------+-------------+-----------+-------------+-------------+---------------+--------------+-----------------------------|
//! | \fbox{h2o} | h2o-w       | Water ice               | cryst       |      0.92 |        0.04 |         2e6 | [[https://ui.adsabs.harvard.edu/abs/2008JGRD..11314220W][Warren+08]]     |              | [[file:lnk_data/h2o-w-Warren2008.lnk][h2o-w-Warren2008.lnk]]        |
//! |            | h2o-a       | Water ice               | amorph      |      0.92 |        0.04 |         2e6 | [[https://ui.adsabs.harvard.edu/abs/1993ApJS...86..713H][Hudgins+93]]    | +Warren      | [[file:lnk_data/h2o-a-Hudgins1993.lnk][h2o-a-Hudgins1993.lnk]]       |
//! | co2        | co2-w       | CO_2 ice                | cryst       |       1.6 |        0.05 |         2e5 | [[https://ui.adsabs.harvard.edu/abs/1986ApOpt..25.2650W][Warren+86]]     | interpolated | [[file:lnk_data/co2-ice-Warren1986.lnk][co2-ice-Warren2008.lnk]]      |
//! | nh3        | nh3-m       | NH_3 ice                | cryst       |      0.75 |        0.14 |         200 | [[https://ui.adsabs.harvard.edu/abs/1984ApOpt..23..541M][Martonchik+83]] | \rho?        | [[file:lnk_data/nh3-m-Martonchik1983.lnk][nh3-m-Martonchik1983.lnk]]    |
//! | co         | co-a        | CO ice                  | amorph      |      0.81 |       *3.8* |       *5.8* | [[https://ui.adsabs.harvard.edu/abs/2006PCCP....8..279P][Palumbo+06]]    |              | [[file:lnk_data/co-a-Palumbo2006.lnk][co-a-Palumbo2006.lnk]]        |
//! |            | co2-a / c   | CO_2 ice                | am / cr     |       1.2 |       *2.5* |        *20* | [[https://ui.adsabs.harvard.edu/abs/2020ApJ...901...52G][Gerakines+20]]  |              | [[file:lnk_data/co2-a-Gerakines2020.lnk][amorph]]/[[file:lnk_data/co2-c-Gerakines2020.lnk][cryst]]                |
//! |            | ch4-a / c   | CH_4 ice                | am / cr     |      0.47 |       *2.0* |        *20* | [[https://ui.adsabs.harvard.edu/abs/2020ApJ...901...52G][Gerakines+20]]  |              | [[file:lnk_data/ch4-a-Gerakines2020.lnk][amorph]]/[[file:lnk_data/ch4-c-Gerakines2020.lnk][cryst]]                |
//! |            | ch3oh-a / c | CH_{3}OH ice            | am / cr     | 0.78/1.02 |       *2.0* |        *24* | [[https://ui.adsabs.harvard.edu/abs/2020ApJ...901...52G][Gerakines+20]]  |              | [[file:lnk_data/ch3oh-a-Gerakines2020.lnk][amorph]]/[[file:lnk_data/ch3oh-c-Gerakines2020.lnk][cryst]]                |

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
