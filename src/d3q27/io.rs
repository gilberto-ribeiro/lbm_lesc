use super::{Lattice, Residuals, ShallowLattice, Simulation};
use crate::io::WriteDataMode;
use colored::*;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::Path;
use std::process;

impl Simulation {
    pub fn build_case_setup() -> io::Result<Simulation> {
        if let Err(e) = crate::io::create_case_directories() {
            eprintln!("Error while creating the case directories: {e}.");
            process::exit(1);
        };
        let case_setup_path =
            Path::new(crate::io::PRE_PROCESSING_PATH).join(crate::io::CASE_SETUP_FILE);
        let case_setup_path_str = case_setup_path.to_str().unwrap();
        let simulation: Simulation;
        if case_setup_path.exists() {
            println!(
                "Reading the case setup file: {}.\n",
                case_setup_path_str.yellow().bold()
            );
            let parameters = crate::io::read_case_setup()?;
            simulation = Simulation::from_setup(parameters);
        } else {
            simulation = Simulation::new();
        }
        if let Err(e) = simulation.create_script_for_residuals_graph() {
            eprintln!("Error while creating script for residuals graph file: {e}.");
            process::exit(1);
        };
        if let Err(e) = simulation.create_script_for_live_residuals_graph() {
            eprintln!("Error while creating script for live residuals graph file: {e}.");
            process::exit(1);
        };
        Ok(simulation)
    }
}

impl Simulation {
    pub fn print_residuals(&self, residuals: &Residuals) {
        if self.time_step % 100 == 0 {
            let duration = self.simulation_time.elapsed().as_secs_f64();
            println!("\n{} {:.2} s.", "Elapsed time:".cyan().bold(), duration);
            println!(
                "\n{:>8} {:>16} {:>16} {:>16} {:>16}\n",
                "step".cyan().bold(),
                "density".cyan().bold(),
                "velocity_x".cyan().bold(),
                "velocity_y".cyan().bold(),
                "velocity_z".cyan().bold()
            );
        }
        println!(
            "{:>8} {:>16.8e} {:>16.8e} {:>16.8e} {:>16.8e}",
            self.time_step,
            residuals.density,
            residuals.velocity[0],
            residuals.velocity[1],
            residuals.velocity[2]
        );
    }

    pub fn write_residuals(&self, residuals: &Residuals) -> io::Result<()> {
        let data_path = Path::new(crate::io::DATA_PATH);
        let path = data_path.join(crate::io::RESIDUALS_FILE);
        let mut file = OpenOptions::new().create(true).append(true).open(path)?;
        if self.time_step == 0 {
            writeln!(
                file,
                "{:>8} {:>16} {:>16} {:>16} {:>16}",
                "step", "density", "velocity_x", "velocity_y", "velocity_z"
            )?;
        }
        writeln!(
            file,
            "{:>8} {:>16.8e} {:>16.8e} {:>16.8e} {:>16.8e}",
            self.time_step,
            residuals.density,
            residuals.velocity[0],
            residuals.velocity[1],
            residuals.velocity[2],
        )?;
        Ok(())
    }

    pub fn write_data(&self, lattice: &ShallowLattice) {
        match &self.write_data_mode {
            WriteDataMode::Frequency(n) => {
                if self.time_step % n == 0 || self.time_step == 0 {
                    println!();
                    self.write_data_from_steps(lattice);
                }
            }
            WriteDataMode::ListOfSteps(list) => {
                if list.contains(&self.time_step) || self.time_step == 0 {
                    println!();
                    self.write_data_from_steps(lattice);
                }
            }
        }
    }

    pub fn write_post_processing_from_each_n_steps<F>(
        &self,
        lattice: &Lattice,
        n: usize,
        function: F,
        file_name: &str,
    ) -> io::Result<()>
    where
        F: Fn(&Lattice) -> Vec<crate::post::PostResult>,
    {
        if self.time_step % n == 0 {
            let post_results = &function(lattice);
            let post_processing_path = Path::new(crate::io::POST_PROCESSING_PATH);
            let path = post_processing_path.join(file_name);
            let mut file = OpenOptions::new().create(true).append(true).open(path)?;
            if self.time_step == 0 {
                write!(file, "{:>8}", "step")?;
                for post_result in post_results {
                    write!(file, " {:>16}", post_result.name)?;
                }
                writeln!(file)?;
            }
            write!(file, "{:>8}", self.time_step)?;
            for post_result in post_results {
                write!(file, " {:>16.8e}", post_result.value)?;
            }
            writeln!(file)?;
        }
        Ok(())
    }

    fn create_script_for_residuals_graph(&self) -> io::Result<()> {
        let post_processing_path = Path::new(crate::io::POST_PROCESSING_PATH);
        let path = post_processing_path.join(crate::io::RESIDUALS_GRAPH_FILE);
        let mut file = File::create(&path)?;
        let path_str = path.to_str().unwrap();
        println!(
            "Creating the residuals graph script file: {}.\n",
            path_str.yellow().bold()
        );
        writeln!(
            file,
            r#"set title "{case_name}"
    set ylabel "Residuals"
    set xlabel "Iterations"
    set grid
    set logscale y
    set yrange [{min_tolerance}:]
    set ytics format "%L"
    set mxtics 5
    set terminal push
    set terminal pngcairo font "courier"
    set output "fig_{case_name_prefix}_residuals.png"
    plot "../data/residuals.dat" u 1:2 t "density" w l,\
    "" u 1:3 t "velocity (x)" w l,\
    "" u 1:4 t "velocity (y)" w l,\
    "" u 1:5 t "velocity (z)" w l
    set terminal pdfcairo font "courier"
    set output "fig_{case_name_prefix}_residuals.pdf"
    replot
    set terminal pop
    set output"#,
            case_name = self.case_name,
            min_tolerance = self.tolerance_density,
            case_name_prefix = self.case_name.replace(" ", "_").to_lowercase(),
        )?;
        Ok(())
    }

    fn create_script_for_live_residuals_graph(&self) -> io::Result<()> {
        let post_processing_path = Path::new(crate::io::POST_PROCESSING_PATH);
        let path = post_processing_path.join(crate::io::LIVE_RESIDUALS_GRAPH_FILE);
        let mut file = File::create(&path)?;
        let path_str = path.to_str().unwrap();
        println!(
            "Creating the live residuals graph script: {}.\n",
            path_str.yellow().bold()
        );
        writeln!(
            file,
            r#"set title "{case_name}"
    set ylabel "Residuals"
    set xlabel "Iterations"
    set grid
    set logscale y
    set yrange [{min_tolerance}:]
    set ytics format "%L"
    set mxtics 5
    set terminal push
    set terminal qt font "courier,12"
    bind "q" "true=0"
    true = 1
    while (true) {{
    plot "../data/residuals.dat" u 1:2 t "density" w l lw 2,\
    "" u 1:3 t "velocity (x)" w l lw 2,\
    "" u 1:4 t "velocity (y)" w l lw 2,\
    "" u 1:5 t "velocity (z)" w l lw 2
    pause 1
    }}
    set terminal pop"#,
            case_name = self.case_name,
            min_tolerance = self.tolerance_density,
        )?;
        Ok(())
    }

    pub fn write_data_from_steps(&self, lattice: &ShallowLattice) {
        let data_path = Path::new(crate::io::DATA_PATH);
        let step_path = data_path.join(&self.time_step.to_string());
        if let Err(e) = fs::create_dir_all(&step_path) {
            eprintln!(
                "Error while creating the the step {} directory: {e}.",
                &self.time_step
            );
            process::exit(1);
        };
        let path = step_path.join(crate::io::DENSITY_FILE);
        println!(
            "\nWriting {} for time step {}.\n",
            crate::io::DENSITY_FILE.yellow().bold(),
            &self.time_step.to_string().yellow().bold()
        );
        if let Err(e) = write_density(lattice, path) {
            eprintln!("Error while writing the density file: {e}.");
            process::exit(1);
        };
        let path = step_path.join(crate::io::VELOCITY_FILE);
        println!(
            "\nWriting {} for time step {}.\n",
            crate::io::VELOCITY_FILE.yellow().bold(),
            &self.time_step.to_string().yellow().bold(),
        );
        if let Err(e) = write_velocity(lattice, path) {
            eprintln!("Error while writing the velocity file: {e}.");
            process::exit(1);
        };
    }

    pub fn write_vtk_from_steps(&self, lattice: &ShallowLattice) {
        let vtk_files_path = Path::new(crate::io::VTK_PATH);
        let case_name = self.case_name.replace(" ", "_").to_lowercase();
        let path_str = format!("{}_results_{:08}.vtk", case_name, &self.time_step);
        let path = vtk_files_path.join(&path_str);
        println!(
            "\nWriting {} for time step {}.\n",
            path_str.yellow().bold(),
            &self.time_step.to_string().yellow().bold()
        );
        if let Err(e) = write_vtk(lattice, path) {
            eprintln!("Error while writing the vtk file: {e}.");
            process::exit(1);
        };
    }

    pub fn write_node_type_vtk(&self, lattice: &Lattice) -> io::Result<()> {
        if lattice.bounce_back_map {
            let vtk_files_path = Path::new(crate::io::VTK_PATH);
            let case_name = self.case_name.replace(" ", "_").to_lowercase();
            let path_str = format!("{}_node_type.vtk", case_name);
            let path = vtk_files_path.join(&path_str);
            println!(
                "Writing node type vtk file: {}.\n",
                path.to_str().unwrap().yellow().bold()
            );
            let mut file = File::create(path)?;
            writeln!(file, "# vtk DataFile Version 3.0")?;
            writeln!(file, "LBM simulation data")?;
            writeln!(file, "ASCII")?;
            writeln!(file, "DATASET STRUCTURED_GRID")?;
            writeln!(
                file,
                "DIMENSIONS {} {} {}",
                lattice.nx, lattice.ny, lattice.nz
            )?;
            writeln!(
                file,
                "POINTS {} float",
                lattice.nx * lattice.ny * lattice.nz
            )?;
            for k in 0..lattice.nz {
                for j in 0..lattice.ny {
                    for i in 0..lattice.nx {
                        writeln!(
                            file,
                            "{:>.4e} {:>.4e} {:>.4e}",
                            lattice.get_node(&[i, j, k]).coordinates[0],
                            lattice.get_node(&[i, j, k]).coordinates[1],
                            lattice.get_node(&[i, j, k]).coordinates[2]
                        )?;
                    }
                }
            }
            writeln!(file, "POINT_DATA {}", lattice.nx * lattice.ny * lattice.nz)?;
            writeln!(file, "SCALARS node_type float 1")?;
            writeln!(file, "LOOKUP_TABLE default")?;
            for k in 0..lattice.nz {
                for j in 0..lattice.ny {
                    for i in 0..lattice.nx {
                        match lattice.get_node(&[i, j, k]).node_type {
                            super::NodeType::Fluid => writeln!(file, "0")?,
                            super::NodeType::Solid => writeln!(file, "1")?,
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

impl Lattice {
    fn write_coordinates(&self) -> io::Result<()> {
        let data_path = Path::new(crate::io::DATA_PATH);
        let path = data_path.join(crate::io::COORDINATES_FILE);
        println!("Writing {}.\n", crate::io::COORDINATES_FILE.yellow().bold());
        let mut file = File::create(path)?;
        writeln!(
            file,
            "{:>8} {:>8} {:>8} {:>16} {:>16} {:>16}",
            "i", "j", "k", "x", "y", "z"
        )?;
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    writeln!(
                        file,
                        "{i:>8} {j:>8} {k:>8} {x:>16.8e} {y:>16.8e} {z:>16.8e}",
                        x = self.get_node(&[i, j, k]).coordinates[0],
                        y = self.get_node(&[i, j, k]).coordinates[1],
                        z = self.get_node(&[i, j, k]).coordinates[2],
                    )?;
                }
            }
        }
        Ok(())
    }

    fn write_case_parameters(&self) -> io::Result<()> {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let case_parameters = self.case_parameters.as_ref().unwrap();
        let length_conversion_factor = case_parameters.length_conversion_factor;
        let time_conversion_factor = case_parameters.time_conversion_factor;
        let density_conversion_factor = case_parameters.density_conversion_factor;
        let velocity_conversion_factor = case_parameters.velocity_conversion_factor;
        let pressure_conversion_factor = case_parameters.pressure_conversion_factor;
        let viscosity_conversion_factor = case_parameters.viscosity_conversion_factor;
        let physical_density = case_parameters.physical_density;
        let reference_pressure = case_parameters.reference_pressure;
        let viscosity = case_parameters.viscosity;
        let physical_viscosity = case_parameters.physical_viscosity;
        let case_parameters_str = format!(
            r#"nx                               = {nx}
ny                               = {ny}
nz                               = {nz}

length_conversion_factor         = {length_conversion_factor:.8e}
time_conversion_factor           = {time_conversion_factor:.8e}
density_conversion_factor        = {density_conversion_factor:.8e}
velocity_conversion_factor       = {velocity_conversion_factor:.8e}
pressure_conversion_factor       = {pressure_conversion_factor:.8e}
viscosity_conversion_factor      = {viscosity_conversion_factor:.8e}

physical_density                 = {physical_density:.8e}
reference_pressure               = {reference_pressure:.8e}

viscosity                        = {viscosity:.8e}
physical_viscosity               = {physical_viscosity:.8e}"#
        );
        let post_processing_path = Path::new(crate::io::POST_PROCESSING_PATH);
        let path = post_processing_path.join(crate::io::CASE_PARAMETERS_FILE);
        println!(
            "Writing the case parameters file: {}.\n",
            path.to_str().unwrap().yellow().bold()
        );
        let mut file = File::create(path)?;
        writeln!(file, "{}", case_parameters_str)?;
        Ok(())
    }
}

impl Lattice {
    pub fn build_case_conditions() -> io::Result<Lattice> {
        let case_conditions_path =
            Path::new(crate::io::PRE_PROCESSING_PATH).join(crate::io::CASE_CONDITIONS_FILE);
        let case_conditions_path_str = case_conditions_path.to_str().unwrap();
        if case_conditions_path.exists() {
            println!(
                "Reading the case conditions file: {}.\n",
                case_conditions_path_str.yellow().bold()
            );
        } else {
            let default_case_conditions = String::from(
                r#"nx                               = 400
ny                               = 100
nz                               = 100

bounce_back_map                  = false

initial_velocity                 = 0.0 0.0 0.0
initial_density                  = 1.0

west_boundary_condition          = dirichlet 1.0 0.1 0.0 0.0
east_boundary_condition          = pressure_outlet 1.0
south_boundary_condition         = no_slip
north_boundary_condition         = no_slip
bottom_boundary_condition        = no_slip
top_boundary_condition           = no_slip
"#,
            );
            let mut file = File::create(&case_conditions_path)?;
            println!(
                "Creating the default case conditions file: {}.\n",
                case_conditions_path_str.yellow().bold()
            );
            write!(file, "{}", default_case_conditions)?;
        }
        let conditions = crate::io::read_case_conditions()?;
        let lattice = Lattice::initialization(conditions);
        lattice.write_coordinates()?;
        lattice.write_case_parameters()?;
        Ok(lattice)
    }
}

pub fn read_bounce_back_map() -> io::Result<Vec<Vec<Vec<i32>>>> {
    let pre_processing_path = Path::new(crate::io::PRE_PROCESSING_PATH);
    let path = pre_processing_path.join(crate::io::BOUNCE_BACK_MAP_FILE);
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    let mut map = Vec::new();
    let mut slice = Vec::new();
    for line in contents.lines() {
        let line = line.trim();
        if line.is_empty() {
            if !slice.is_empty() {
                map.push(slice);
                slice = Vec::new();
            }
            continue;
        }
        let mut row = Vec::new();
        for value in line.split_whitespace() {
            row.push(value.parse::<i32>().unwrap());
        }
        slice.push(row);
    }
    if !slice.is_empty() {
        map.push(slice);
    }
    map = y_axis_reverse_matrix(map);
    Ok(map)
}

fn write_density<P>(lattice: &ShallowLattice, path: P) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut file = File::create(path)?;
    writeln!(file, "{:>16}", "density")?;
    for k in 0..lattice.nz {
        for j in 0..lattice.ny {
            for i in 0..lattice.nx {
                writeln!(
                    file,
                    "{density:>16.8e}",
                    density = lattice.get_node(&[i, j, k]).density,
                )?;
            }
        }
    }
    Ok(())
}

fn write_velocity<P>(lattice: &ShallowLattice, path: P) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut file = File::create(path)?;
    writeln!(
        file,
        "{:>16} {:>16} {:>16}",
        "velocity_x", "velocity_y", "velocity_z"
    )?;
    for k in 0..lattice.nz {
        for j in 0..lattice.ny {
            for i in 0..lattice.nx {
                writeln!(
                    file,
                    "{velocity_x:>16.8e} {velocity_y:>16.8e} {velocity_z:>16.8e}",
                    velocity_x = lattice.get_node(&[i, j, k]).velocity[0],
                    velocity_y = lattice.get_node(&[i, j, k]).velocity[1],
                    velocity_z = lattice.get_node(&[i, j, k]).velocity[2],
                )?;
            }
        }
    }
    Ok(())
}

fn write_vtk<P>(lattice: &ShallowLattice, path: P) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut file = File::create(path)?;
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "LBM simulation data")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET STRUCTURED_GRID")?;
    writeln!(
        file,
        "DIMENSIONS {} {} {}",
        lattice.nx, lattice.ny, lattice.nz
    )?;
    writeln!(
        file,
        "POINTS {} float",
        lattice.nx * lattice.ny * lattice.nz
    )?;
    for k in 0..lattice.nz {
        for j in 0..lattice.ny {
            for i in 0..lattice.nx {
                writeln!(
                    file,
                    "{:>.4e} {:>.4e} {:>.4e}",
                    lattice.get_node(&[i, j, k]).coordinates[0],
                    lattice.get_node(&[i, j, k]).coordinates[1],
                    lattice.get_node(&[i, j, k]).coordinates[2]
                )?;
            }
        }
    }
    writeln!(file, "POINT_DATA {}", lattice.nx * lattice.ny * lattice.nz)?;
    writeln!(file, "SCALARS density float 1")?;
    writeln!(file, "LOOKUP_TABLE default")?;
    for k in 0..lattice.nz {
        for j in 0..lattice.ny {
            for i in 0..lattice.nx {
                writeln!(file, "{:>.8e}", lattice.get_node(&[i, j, k]).density)?;
            }
        }
    }
    writeln!(file, "VECTORS velocity float")?;
    for k in 0..lattice.nz {
        for j in 0..lattice.ny {
            for i in 0..lattice.nx {
                writeln!(
                    file,
                    "{:>.8e} {:>.8e} {:>.8e}",
                    lattice.get_node(&[i, j, k]).velocity[0],
                    lattice.get_node(&[i, j, k]).velocity[1],
                    lattice.get_node(&[i, j, k]).velocity[2]
                )?;
            }
        }
    }
    Ok(())
}

fn y_axis_reverse_matrix(matrix: Vec<Vec<Vec<i32>>>) -> Vec<Vec<Vec<i32>>> {
    let mut reversed_matrix = matrix;
    for slice in reversed_matrix.iter_mut() {
        slice.reverse();
    }
    reversed_matrix
}
