use anyhow::Result;

use kappa::cli;
use kappa::opac::KappaError;
use kappa::opac::run;

fn main() -> Result<(), KappaError> {
    let mut kpc = cli::launch()?;
    run(&mut kpc)?;

    Ok(())
}
