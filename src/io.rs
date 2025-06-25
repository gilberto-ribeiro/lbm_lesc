use colored::*;
use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::Path;
use std::process;
use std::time::Duration;

pub const DATA_PATH: &'static str = "./data";

pub const PRE_PROCESSING_PATH: &'static str = "./pre_processing";

pub const CASE_SETUP_FILE: &'static str = "case_setup.jou";

pub const CASE_CONDITIONS_FILE: &'static str = "case_conditions.jou";

pub const CASE_PARAMETERS_FILE: &'static str = "case_parameters.jou";

pub const POST_PROCESSING_PATH: &'static str = "./post_processing";

pub const VTK_PATH: &'static str = "./post_processing/vtk_files";

pub const COORDINATES_FILE: &'static str = "coordinates.dat";

pub const DENSITY_FILE: &'static str = "density.dat";

pub const VELOCITY_FILE: &'static str = "velocity.dat";

pub const RESIDUALS_FILE: &'static str = "residuals.dat";

pub const BOUNCE_BACK_MAP_FILE: &'static str = "map.dat";

pub const RESIDUALS_GRAPH_FILE: &'static str = "gr_residuals.gp";

pub const LIVE_RESIDUALS_GRAPH_FILE: &'static str = "live_residuals.gp";

#[derive(Clone)]
pub enum WriteDataMode {
    Frequency(usize),

    ListOfSteps(Vec<usize>),
}

pub fn create_case_directories() -> io::Result<()> {
    let list_of_paths = [
        DATA_PATH,
        PRE_PROCESSING_PATH,
        POST_PROCESSING_PATH,
        VTK_PATH,
    ];
    for path_str in list_of_paths {
        let path = Path::new(path_str);
        if !path.exists() {
            println!("Creating the {} path.\n", path_str.yellow().bold());
            if let Err(e) = fs::create_dir(path) {
                eprintln!("Error while creating the {path_str} path: {e}.");
                process::exit(1);
            };
        } else {
            println!("The {} path already exists.\n", path_str.yellow().bold());
        }
    }
    Ok(())
}

pub fn read_case_setup() -> io::Result<HashMap<String, String>> {
    let path = Path::new(PRE_PROCESSING_PATH).join(CASE_SETUP_FILE);
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(extract_parameters(contents))
}

pub fn read_case_conditions() -> io::Result<HashMap<String, String>> {
    let path = Path::new(PRE_PROCESSING_PATH).join(CASE_CONDITIONS_FILE);
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(extract_parameters(contents))
}

pub fn read_case_parameters() -> io::Result<HashMap<String, String>> {
    let path = Path::new(POST_PROCESSING_PATH).join(CASE_PARAMETERS_FILE);
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(extract_parameters(contents))
}

fn extract_parameters(contents: String) -> HashMap<String, String> {
    contents
        .lines()
        .map(|line| line.trim())
        .filter(|line| !line.starts_with("#"))
        .filter(|line| !line.is_empty())
        .map(|line| {
            let mut parts = line.splitn(2, "=");
            let key = parts.next().unwrap().trim().to_string();
            let value = parts.next().unwrap_or("").trim().to_string();
            (key, value)
        })
        .collect::<HashMap<String, String>>()
}

pub fn write_inside_loop_elapsed_time(
    elapsed_times: &[(&str, Duration)],
    time_step: &usize,
) -> io::Result<()> {
    let path = Path::new(POST_PROCESSING_PATH).join("benchmark_elapsed_time.dat");
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    if *time_step == 0 {
        write!(file, "{:>8}", "step")?;
        for (key, _) in elapsed_times {
            write!(file, " {:>16}", key)?;
        }
        writeln!(file)?;
    }
    write!(file, "{:>8}", time_step)?;
    for (_, value) in elapsed_times {
        write!(file, " {:>16.8e}", value.as_secs_f64())?;
    }
    writeln!(file)?;
    Ok(())
}
