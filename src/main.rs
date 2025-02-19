use anyhow::Result;

use kappa::cli;

fn main() -> Result<()> {
    cli::run()?;

    Ok(())
}
