use anyhow::Result;

use kappa::cli;
use kappa::opac::run;
use kappa::opac::KappaError;

fn main() -> Result<(), KappaError> {
    let mut kpc = cli::launch()?;
    run(&mut kpc)?;

    Ok(())
}
