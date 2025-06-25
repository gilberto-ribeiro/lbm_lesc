use crate::d2q9::D;
use crate::d2q9::{Lattice, ShallowLattice};
use crate::global_variables::*;
use colored::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

pub struct CaseParameters {
    pub nx: usize,
    pub ny: usize,
    pub length_conversion_factor: Float,
    pub time_conversion_factor: Float,
    pub density_conversion_factor: Float,
    pub velocity_conversion_factor: Float,
    pub pressure_conversion_factor: Float,
    pub viscosity_conversion_factor: Float,
    pub physical_density: Float,
    pub reference_pressure: Float,
    pub viscosity: Float,
    pub physical_viscosity: Float,
}

impl CaseParameters {
    fn from_file() -> io::Result<CaseParameters> {
        let parameters = crate::io::read_case_parameters().unwrap();
        Ok(CaseParameters {
            nx: parameters["nx"].parse().unwrap(),
            ny: parameters["ny"].parse().unwrap(),
            length_conversion_factor: parameters["length_conversion_factor"].parse().unwrap(),
            time_conversion_factor: parameters["time_conversion_factor"].parse().unwrap(),
            density_conversion_factor: parameters["density_conversion_factor"].parse().unwrap(),
            velocity_conversion_factor: parameters["velocity_conversion_factor"].parse().unwrap(),
            pressure_conversion_factor: parameters["pressure_conversion_factor"].parse().unwrap(),
            viscosity_conversion_factor: parameters["viscosity_conversion_factor"].parse().unwrap(),
            physical_density: parameters["physical_density"].parse().unwrap(),
            reference_pressure: parameters["reference_pressure"].parse().unwrap(),
            viscosity: parameters["viscosity"].parse().unwrap(),
            physical_viscosity: parameters["physical_viscosity"].parse().unwrap(),
        })
    }

    pub fn new(lattice: &Lattice) -> CaseParameters {
        let nx = lattice.nx;
        let ny = lattice.ny;
        let length_conversion_factor = lattice.delta_x;
        let time_conversion_factor = lattice.delta_t;
        let density_conversion_factor = lattice.physical_density;
        let velocity_conversion_factor = lattice.delta_x / lattice.delta_t;
        let pressure_conversion_factor =
            density_conversion_factor * velocity_conversion_factor * velocity_conversion_factor;
        let viscosity_conversion_factor = lattice.delta_x * lattice.delta_x / lattice.delta_t;
        let physical_density = lattice.physical_density;
        let reference_pressure = lattice.reference_pressure;
        let viscosity = CS_2 * (lattice.tau - 0.5);
        let physical_viscosity = viscosity_conversion_factor * viscosity;
        CaseParameters {
            nx,
            ny,
            length_conversion_factor,
            time_conversion_factor,
            density_conversion_factor,
            velocity_conversion_factor,
            pressure_conversion_factor,
            viscosity_conversion_factor,
            physical_density,
            reference_pressure,
            viscosity,
            physical_viscosity,
        }
    }
}

impl ShallowLattice {
    pub fn read_coordinates(&mut self) {
        let data_path = Path::new(crate::io::DATA_PATH);
        let path = data_path.join(crate::io::COORDINATES_FILE);
        let file = File::open(path).expect("Unable to open file");
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        lines.next();
        lines.zip(self.nodes.iter_mut()).for_each(|(line, node)| {
            let line = line.expect("Unable to read line");
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 2 * D {
                let coordinates: [Float; D] = [
                    parts[2].parse().expect("Invalid coordinates value"),
                    parts[3].parse().expect("Invalid coordinates value"),
                ];
                node.coordinates = coordinates;
            }
        });
    }

    pub fn read_density(&mut self, time_step: usize) {
        let data_path = Path::new(crate::io::DATA_PATH);
        let step_path = data_path.join(time_step.to_string());
        let path = step_path.join(crate::io::DENSITY_FILE);
        let file = File::open(path).expect("Unable to open file");
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        lines.next();
        lines.zip(self.nodes.iter_mut()).for_each(|(line, node)| {
            let line = line.expect("Unable to read line");
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 1 {
                let density: Float = parts[0].parse().expect("Invalid density value");
                node.density = density;
            }
        });
    }

    pub fn read_velocity(&mut self, time_step: usize) {
        let data_path = Path::new(crate::io::DATA_PATH);
        let step_path = data_path.join(time_step.to_string());
        let path = step_path.join(crate::io::VELOCITY_FILE);
        let file = File::open(path).expect("Unable to open file");
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        lines.next();
        lines.zip(self.nodes.iter_mut()).for_each(|(line, node)| {
            let line = line.expect("Unable to read line");
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == D {
                let velocity: [Float; D] = [
                    parts[0].parse().expect("Invalid velocity value"),
                    parts[1].parse().expect("Invalid velocity value"),
                ];
                node.velocity = velocity;
            }
        });
    }
}

impl ShallowLattice {
    pub fn from_data(time_step: usize, case_parameters: &CaseParameters) -> ShallowLattice {
        let mut lattice =
            ShallowLattice::new(case_parameters.nx, case_parameters.ny, 1.0, [0.0; D]);
        lattice.read_coordinates();
        lattice.read_density(time_step);
        lattice.read_velocity(time_step);
        lattice
    }
}

fn write_physical_vtk<P>(
    lattice: &ShallowLattice,
    path: P,
    case_parameters: &CaseParameters,
) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let point_data = lattice.nx * lattice.ny;
    let mut file = File::create(path)?;
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "LBM simulation data")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET STRUCTURED_GRID")?;
    writeln!(file, "DIMENSIONS {} {} 1", lattice.nx, lattice.ny)?;
    writeln!(file, "POINTS {} float", point_data)?;
    lattice.nodes.iter().for_each(|node| {
        writeln!(
            file,
            "{:>.6e} {:>.6e} 0.0",
            node.coordinates[0], node.coordinates[1]
        )
        .unwrap();
    });
    writeln!(file, "POINT_DATA {}", point_data)?;
    writeln!(file, "SCALARS lattice_density float 1")?;
    writeln!(file, "LOOKUP_TABLE default")?;
    lattice.nodes.iter().for_each(|node| {
        let density = node.density;
        writeln!(file, "{density:>.6e}").unwrap();
    });
    writeln!(file, "VECTORS lattice_velocity float")?;
    lattice.nodes.iter().for_each(|node| {
        let velocity_x = node.velocity[0];
        let velocity_y = node.velocity[1];
        writeln!(file, "{velocity_x:>.6e} {velocity_y:>.6e} 0.0").unwrap();
    });
    writeln!(file, "SCALARS physical_pressure float 1")?;
    writeln!(file, "LOOKUP_TABLE default")?;
    lattice.nodes.iter().for_each(|node| {
        let density_prime = node.density - LATTICE_DENSITY;
        let pressure_prime = CS_2 * density_prime;
        let physical_pressure = case_parameters.reference_pressure
            + pressure_prime * case_parameters.pressure_conversion_factor;
        writeln!(file, "{physical_pressure:>.6e}").unwrap();
    });
    writeln!(file, "VECTORS physical_velocity float")?;
    lattice.nodes.iter().for_each(|node| {
        let physical_velocity_x = node.velocity[0] * case_parameters.velocity_conversion_factor;
        let physical_velocity_y = node.velocity[1] * case_parameters.velocity_conversion_factor;
        writeln!(
            file,
            "{physical_velocity_x:>.6e} {physical_velocity_y:>.6e} 0.0"
        )
        .unwrap();
    });
    Ok(())
}

fn read_data_directory() -> io::Result<Vec<usize>> {
    let mut time_steps = Vec::new();
    let path = Path::new(crate::io::DATA_PATH);
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        if entry.file_type()?.is_dir() {
            if let Ok(time_step) = entry.file_name().to_string_lossy().parse::<usize>() {
                time_steps.push(time_step);
            }
        }
    }
    time_steps.sort_unstable();
    Ok(time_steps)
}

pub fn run_vtk_post_processing() -> io::Result<()> {
    let case_setup = crate::io::read_case_setup()?;
    let case_name = case_setup["case_name"]
        .clone()
        .replace(" ", "_")
        .to_lowercase();
    let case_parameters = CaseParameters::from_file()?;
    let time_steps = read_data_directory()?;
    println!("{:?}", time_steps);
    time_steps.par_iter().for_each(|&time_step| {
        let lattice = ShallowLattice::from_data(time_step, &case_parameters);
        let path_str = format!("{case_name}_{:08}.vtk", time_step);
        let path = Path::new(crate::io::VTK_PATH).join(&path_str);
        println!(
            "Writing {} for time step {}.\n",
            path_str.yellow().bold(),
            time_step.to_string().yellow().bold()
        );
        if let Err(e) = write_physical_vtk(&lattice, path, &case_parameters) {
            eprintln!("Error while writing the vtk file: {e}.");
            std::process::exit(1);
        };
    });
    Ok(())
}
