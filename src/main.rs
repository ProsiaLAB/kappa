use anyhow::Result;

use kappa::{cli, opac::KappaError};

fn main() -> Result<(), KappaError> {
    cli::run()?;

    Ok(())
}
