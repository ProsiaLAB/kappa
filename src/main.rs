use anyhow::Result;

use kappa::{cli, opac::KappaError};

fn main() -> Result<(), KappaError> {
    cli::launch()?;

    Ok(())
}
