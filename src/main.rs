use anyhow::Result;

use kappa::cli;
use kappa::opac::KappaError;

fn main() -> Result<(), KappaError> {
    let mut kpc = cli::launch()?;
    kpc.prepare_inputs()?;
    kpc.initialize()?;
    kpc.run()?;

    Ok(())
}
