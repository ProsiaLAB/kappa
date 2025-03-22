//! Read the `lnk` files

use anyhow::Result;

use std::io::BufRead;

use crate::opac::Component;

pub fn read_lnk_file(file: &str, rho_in: Option<f64>) -> Result<Component> {
    let file = std::fs::File::open(file)?;
    let reader = std::io::BufReader::new(file);
    let all_lines: Vec<String> = reader.lines().map_while(Result::ok).collect();
    let mut lines_iter = all_lines.iter();

    let mut name = String::new();
    let mut class = String::new();
    let mut state = String::new();

    let header = lines_iter
        .find_map(|line| {
            let trimmed = line.trim();
            if trimmed.starts_with('#') {
                if let Some((key, value)) = trimmed.trim_start_matches('#').trim().split_once(':') {
                    match key.trim() {
                        "Name" => name = value.trim().to_string(),
                        "Class" => class = value.trim().to_string(),
                        "State" => state = value.trim().to_string(),
                        _ => {}
                    }
                }
                None
            } else if !trimmed.is_empty() {
                Some(trimmed.to_string())
            } else {
                None
            }
        })
        .expect("Missing header in the size distribution file");

    // Extract header values: first value is nlam, second is rho
    let mut header_parts = header.split_whitespace();
    let size: usize = header_parts
        .next()
        .expect("Missing nlam in the size distribution file")
        .parse()
        .expect("Invalid nlam in the size distribution file");
    let rho: f64 = match header_parts
        .next()
        .and_then(|s| s.parse::<f64>().ok())
        .or(rho_in)
    {
        Some(r) => r,
        None => panic!("Missing rho in the size distribution file"),
    };

    // Read the rest of the file
    let mut l0 = vec![0.0; size];
    let mut n0 = vec![0.0; size];
    let mut k0 = vec![0.0; size];

    let mut iter = lines_iter
        .filter(|line| !line.starts_with('#')) // Ignore remaining comments
        .map(|line| {
            let mut values = line
                .split_whitespace()
                .map(|s| s.parse::<f64>().expect("Invalid number in data"));
            (
                values.next().expect("Missing l"),
                values.next().expect("Missing n"),
                values.next().expect("Missing k"),
            )
        });

    for ((l, n), k) in l0.iter_mut().zip(n0.iter_mut()).zip(k0.iter_mut()) {
        let (l_val, n_val, k_val) = iter.next().expect("Too few rows in the file");
        *l = l_val;
        *n = n_val;
        *k = k_val;
    }

    regrid_data(&mut l0, &mut n0, &mut k0);

    let component = Component {
        name,
        class,
        state,
        rho,
        size,
        l0,
        n0,
        k0,
    };

    Ok(component)
}

fn regrid_data(l0: &mut [f64], n0: &mut [f64], k0: &mut [f64]) {
    // Check if we need to reverse the arrays
    if l0[l0.len() - 1] < l0[0] {
        l0.reverse();
        n0.reverse();
        k0.reverse();
    }

    // Disallow values <= 0 to accomodate the log-log interpolation
    for (n, k) in n0.iter_mut().zip(k0.iter_mut()) {
        if *n <= 0.0 {
            *n = 1e-10;
        }
        if *k <= 0.0 {
            *k = 1e-10;
        }
    }

    let x0 = l0[0];
    let y01 = n0[0];
    let y02 = k0[0];
    let _wp = (1.0 - y01) / x0.powi(2);
    let _gamma = y02 / x0.powi(3);
}

pub fn read_sizedis_file(_file: &str) -> Result<()> {
    // let file = std::fs::File::open(file)?;
    // let reader = std::io::BufReader::new(file);

    // // Skip comments
    // let mut lines_iter = reader.lines().map_while(Result::ok);
    // for line in lines_iter.by_ref() {
    //     if !line.starts_with('#') || !line.starts_with("!") || !line.starts_with("*") {
    //         let na = line
    //             .trim()
    //             .parse::<usize>()
    //             .expect("Invalid number of size bins");

    //     }
    // }

    todo!()

    // Ok(())
}

pub fn read_wavelength_grid(_file: &str) -> Result<()> {
    todo!()
}
