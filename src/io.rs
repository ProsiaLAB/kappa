//! Read the `lnk` files

use anyhow::Result;

use std::io::BufRead;

type LNKData = (Vec<f64>, Vec<f64>, Vec<f64>);

pub fn read_lnk_file(file: &str) -> Result<(LNKData, usize, f64)> {
    let file = std::fs::File::open(file)?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();
    let header = lines
        .find_map(|line| {
            let line = line.ok()?;
            if !line.starts_with('#') {
                Some(line)
            } else {
                None
            }
        })
        .expect("Missing header in the size distribution file");
    // Extract header values: first value is nlam, second is rho
    let mut header_parts = header.split_whitespace();
    let nlam: usize = header_parts
        .next()
        .expect("Missing nlam in the size distribution file")
        .parse()
        .expect("Invalid nlam in the size distribution file");
    let rho: f64 = header_parts
        .next()
        .expect("Missing rho in the size distribution file")
        .parse()
        .expect("Invalid rho in the size distribution file");

    // Read the rest of the file
    let mut l_vec = vec![0.0; nlam];
    let mut n_vec = vec![0.0; nlam];
    let mut k_vec = vec![0.0; nlam];

    let mut l_iter = l_vec.iter_mut();
    let mut n_iter = n_vec.iter_mut();
    let mut k_iter = k_vec.iter_mut();

    for line in lines
        .map_while(Result::ok)
        .filter(|line| !line.starts_with('#'))
    {
        let mut values = line
            .split_whitespace()
            .map(|s| s.parse::<f64>().expect("Invalid number in data"));
        *l_iter.next().expect("Too many rows") = values.next().expect("Missing l");
        *n_iter.next().expect("Too many rows") = values.next().expect("Missing n");
        *k_iter.next().expect("Too many rows") = values.next().expect("Missing k");
    }

    Ok(((l_vec, n_vec, k_vec), nlam, rho))
}
