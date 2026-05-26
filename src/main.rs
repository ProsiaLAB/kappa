use anyhow::Result;

use kappa::cli;
use kappa::opac::KappaError;
use kappa::opac::run;

fn main() -> Result<(), KappaError> {
    let mut kpc = cli::launch()?;
    kpc.prepare_inputs()?;
    kpc.initialize()?;
    run(&kpc)?;

    Ok(())
}
