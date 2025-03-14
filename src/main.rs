use anyhow::Result;

use kappa::{cli, opac::KappaError};

fn main() -> Result<(), KappaError> {
    cli::launch()?;
    println!("No arguments provided. Running in default mode.");

    Ok(())
}
