//! Read the `lnk` files

use anyhow::anyhow;
use anyhow::Result;

use std::io::BufRead;

use crate::opac::Component;
use crate::types::RVector;

pub fn read_lnk_file(file: &str, rho_in: Option<f64>) -> Result<Component> {
    let file = std::fs::File::open(file)?;
    let reader = std::io::BufReader::new(file);
    let all_lines: Vec<String> = reader.lines().map_while(Result::ok).collect();
    let mut lines_iter = all_lines.iter();

    let mut name = String::new();
    let mut class = String::new();
    let mut state = String::new();

    let header: String = lines_iter
        .find_map(|line| {
            let trimmed = line.trim();
            if trimmed.starts_with('#') {
                if let Some((key, value)) = trimmed.trim_start_matches('#').trim().split_once(':') {
                    match key.trim() {
                        "Name" => name = value.trim().into(),
                        "Class" => class = value.trim().into(),
                        "State" => state = value.trim().into(),
                        _ => {}
                    }
                }
                None
            } else if !trimmed.is_empty() {
                Some(trimmed.into())
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

    let component = Component {
        name,
        class,
        state,
        rho,
        size,
        l0: RVector::from_vec(l0),
        n0: RVector::from_vec(n0),
        k0: RVector::from_vec(k0),
    };

    Ok(component)
}

pub fn read_sizedis_file(file: &str) -> Result<(usize, [f64; 3])> {
    let file = std::fs::File::open(file)?;
    let reader = std::io::BufReader::new(file);

    let all_lines: Vec<String> = reader.lines().map_while(Result::ok).collect();
    let mut lines_iter = all_lines.iter();
    let na: usize = lines_iter
        .find_map(|line| {
            let trimmed = line.trim();
            if trimmed.starts_with('#') || trimmed.starts_with('!') || trimmed.starts_with('*') {
                None
            } else if !trimmed.is_empty() {
                // Parse usize value
                Some(
                    trimmed
                        .parse::<usize>()
                        .expect("Invalid header in the size distribution file"),
                )
            } else {
                None
            }
        })
        .expect("Missing header in the size distribution file");

    let mut tot: [f64; 3] = [0.0; 3];
    let mut totn = 0.0;

    for line in lines_iter {
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with('!')
            || trimmed.starts_with('*')
        {
            continue;
        }
        let mut parts = trimmed.split_whitespace();
        if let (Some(a), Some(b)) = (parts.next(), parts.next()) {
            if let (Ok(r), Ok(nr)) = (a.parse::<f64>(), b.parse::<f64>()) {
                totn += nr;
                tot[0] += nr * r;
                tot[1] += nr * r.powi(2);
                tot[2] += nr * r.powi(3);
            }
        }
    }

    let sdmns: [f64; 3] = [
        tot[0] / totn,
        (tot[1] / totn).sqrt(),
        (tot[2] / totn).powf(1.0 / 3.0),
    ];

    Ok((na, sdmns))
}

pub fn read_wavelength_grid(file: &str) -> Result<(usize, RVector)> {
    let file = std::fs::File::open(file)?;
    let reader = std::io::BufReader::new(file);

    let all_lines: Vec<String> = reader.lines().map_while(Result::ok).collect();
    let mut lines_iter = all_lines.iter().skip_while(|line| line.starts_with('#'));

    let nlam = lines_iter
        .next()
        .ok_or_else(|| anyhow!("Missing header in the size distribution file"))?
        .trim()
        .parse::<usize>()?;

    let lam = RVector::from_vec(
        lines_iter
            .map(|line| line.trim().parse::<f64>())
            .collect::<Result<Vec<f64>, _>>()
            .map_err(|e| anyhow!(e))?,
    );

    Ok((nlam, lam))
}
