pub mod vtk;

use super::*;
use crate::global_variables::*;
use crate::post::PostResult;

pub fn compute_mean_velocities(lattice: &Lattice) -> Vec<PostResult> {
    let ux_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.velocity[0])
        .sum::<Float>();
    let uy_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.velocity[1])
        .sum::<Float>();
    let uz_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.velocity[2])
        .sum::<Float>();
    let u_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| {
            (node.velocity[0] * node.velocity[0]
                + node.velocity[1] * node.velocity[1]
                + node.velocity[2] * node.velocity[2])
                .sqrt()
        })
        .sum::<Float>();
    let number_of_fluid_nodes = lattice.fluid_nodes.as_ref().unwrap().len() as Float;
    let ux_mean = ux_sum / number_of_fluid_nodes;
    let uy_mean = uy_sum / number_of_fluid_nodes;
    let uz_mean = uz_sum / number_of_fluid_nodes;
    let u_mean = u_sum / number_of_fluid_nodes;
    let u_result: PostResult = PostResult::new(
        "mean_velocity".to_string(),
        "mean velocity (magnitude)".to_string(),
        u_mean,
        None,
    );
    let ux_result: PostResult = PostResult::new(
        "mean_velocity_x".to_string(),
        "mean velocity (x)".to_string(),
        ux_mean,
        None,
    );
    let uy_result: PostResult = PostResult::new(
        "mean_velocity_y".to_string(),
        "mean velocity (y)".to_string(),
        uy_mean,
        None,
    );
    let uz_result: PostResult = PostResult::new(
        "mean_velocity_z".to_string(),
        "mean velocity (z)".to_string(),
        uz_mean,
        None,
    );
    vec![u_result, ux_result, uy_result, uz_result]
}

pub fn compute_mean_density(lattice: &Lattice) -> Vec<PostResult> {
    let rho_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.density)
        .sum::<Float>();
    let number_of_fluid_nodes = lattice.fluid_nodes.as_ref().unwrap().len() as Float;
    let rho_mean = rho_sum / number_of_fluid_nodes;
    let rho_result: PostResult = PostResult::new(
        "mean_density".to_string(),
        "mean density".to_string(),
        rho_mean,
        None,
    );
    vec![rho_result]
}

pub fn compute_porosity(lattice: &Lattice) -> Vec<PostResult> {
    let number_of_solid_nodes = lattice.solid_nodes.as_ref().unwrap().len() as Float;
    let number_of_fluid_nodes = lattice.fluid_nodes.as_ref().unwrap().len() as Float;
    let porosity = number_of_fluid_nodes / (number_of_fluid_nodes + number_of_solid_nodes);
    let number_of_solid_nodes_result: PostResult = PostResult::new(
        "n_solid_nodes".to_string(),
        "number of solid nodes".to_string(),
        number_of_solid_nodes,
        None,
    );
    let number_of_fluid_nodes_result: PostResult = PostResult::new(
        "n_fluid_nodes".to_string(),
        "number of fluid nodes".to_string(),
        number_of_fluid_nodes,
        None,
    );
    let porosity_result: PostResult = PostResult::new(
        "porosity".to_string(),
        "porosity".to_string(),
        porosity,
        None,
    );
    vec![
        number_of_solid_nodes_result,
        number_of_fluid_nodes_result,
        porosity_result,
    ]
}

pub fn compute_tortuosities(lattice: &Lattice) -> Vec<PostResult> {
    let ux_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.velocity[0])
        .sum::<Float>();
    let uy_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.velocity[1])
        .sum::<Float>();
    let uz_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.velocity[2])
        .sum::<Float>();
    let u_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| {
            (node.velocity[0] * node.velocity[0]
                + node.velocity[1] * node.velocity[1]
                + node.velocity[2] * node.velocity[2])
                .sqrt()
        })
        .sum::<Float>();
    let tortuosity_x = u_sum / ux_sum;
    let tortuosity_y = u_sum / uy_sum;
    let tortuosity_z = u_sum / uz_sum;
    let tortuosity_x_result: PostResult = PostResult::new(
        "tortuosity_x".to_string(),
        "tortuosity (x)".to_string(),
        tortuosity_x,
        None,
    );
    let tortuosity_y_result: PostResult = PostResult::new(
        "tortuosity_y".to_string(),
        "tortuosity (y)".to_string(),
        tortuosity_y,
        None,
    );
    let tortuosity_z_result: PostResult = PostResult::new(
        "tortuosity_z".to_string(),
        "tortuosity (z)".to_string(),
        tortuosity_z,
        None,
    );
    vec![
        tortuosity_x_result,
        tortuosity_y_result,
        tortuosity_z_result,
    ]
}

pub fn compute_max_velocity(lattice: &Lattice) -> Vec<PostResult> {
    let max_velocity = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| {
            (node.velocity[0] * node.velocity[0]
                + node.velocity[1] * node.velocity[1]
                + node.velocity[2] * node.velocity[2])
                .sqrt()
        })
        .reduce_with(|a, b| a.max(b))
        .unwrap_or(0.0);
    let max_velocity_result: PostResult = PostResult::new(
        "max_velocity".to_string(),
        "maximum velocity".to_string(),
        max_velocity,
        None,
    );
    vec![max_velocity_result]
}

pub fn compute_mean_pressures_x(lattice: &Lattice) -> Vec<PostResult> {
    let case_parameters = lattice.case_parameters.as_ref().unwrap();
    let pressure_inlet_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .filter(|node| node.index[0] == 0)
        .map(|node| {
            let density_prime = node.density - LATTICE_DENSITY;
            let pressure_prime = CS_2 * density_prime;
            let physical_pressure = case_parameters.reference_pressure
                + pressure_prime * case_parameters.pressure_conversion_factor;
            physical_pressure
        })
        .sum::<Float>();
    let pressure_outlet_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .filter(|node| node.index[0] == lattice.nx - 1)
        .map(|node| {
            let density_prime = node.density - LATTICE_DENSITY;
            let pressure_prime = CS_2 * density_prime;
            let physical_pressure = case_parameters.reference_pressure
                + pressure_prime * case_parameters.pressure_conversion_factor;
            physical_pressure
        })
        .sum::<Float>();
    let number_of_liquid_nodes_inlet = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .filter(|node| node.index[0] == 0)
        .count() as Float;
    let number_of_liquid_nodes_outlet = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .filter(|node| node.index[0] == lattice.nx - 1)
        .count() as Float;
    let pressure_inlet_mean = pressure_inlet_sum / number_of_liquid_nodes_inlet;
    let pressure_outlet_mean = pressure_outlet_sum / number_of_liquid_nodes_outlet;
    let pressure_inlet_mean_result: PostResult = PostResult::new(
        "pressure_inlet".to_string(),
        "pressure inlet mean".to_string(),
        pressure_inlet_mean,
        None,
    );
    let pressure_outlet_mean_result: PostResult = PostResult::new(
        "pressure_outlet".to_string(),
        "pressure outlet mean".to_string(),
        pressure_outlet_mean,
        None,
    );
    vec![pressure_inlet_mean_result, pressure_outlet_mean_result]
}

pub fn compute_permeability_x(lattice: &Lattice) -> Vec<PostResult> {
    let u_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| {
            (node.velocity[0] * node.velocity[0]
                + node.velocity[1] * node.velocity[1]
                + node.velocity[2] * node.velocity[2])
                .sqrt()
        })
        .sum::<Float>();
    let mut pressure_diff_1 = 0.0;
    let mut pressure_diff_2 = 0.0;
    let mut pressure_diff_3 = 0.0;
    for k in 0..lattice.nz {
        for j in 0..lattice.ny {
            let rho_in_1 = lattice.get_node(&[0, j, k]).density;
            let rho_out_1 = lattice.get_node(&[lattice.nx - 1, j, k]).density;
            let rho_in_2 = lattice.get_node(&[1, j, k]).density;
            let rho_out_2 = lattice.get_node(&[lattice.nx - 2, j, k]).density;
            let pressure_in_1 = CS_2 * rho_in_1;
            let pressure_out_1 = CS_2 * rho_out_1;
            let pressure_in_2 = CS_2 * rho_in_2;
            let pressure_out_2 = CS_2 * rho_out_2;
            pressure_diff_1 += pressure_in_1 - pressure_out_1;
            pressure_diff_2 += pressure_in_1 - pressure_out_2;
            pressure_diff_3 += pressure_in_2 - pressure_out_2;
        }
    }
    let permeability_x_1 = u_sum / pressure_diff_1;
    let permeability_x_2 = u_sum / pressure_diff_2;
    let permeability_x_3 = u_sum / pressure_diff_3;
    let permeability_x_1_result: PostResult = PostResult::new(
        "permeability_x_1".to_string(),
        "permeability_x (1)".to_string(),
        permeability_x_1,
        None,
    );
    let permeability_x_2_result: PostResult = PostResult::new(
        "permeability_x_2".to_string(),
        "permeability_x (2)".to_string(),
        permeability_x_2,
        None,
    );
    let permeability_x_3_result: PostResult = PostResult::new(
        "permeability_x_3".to_string(),
        "permeability_x (3)".to_string(),
        permeability_x_3,
        None,
    );
    vec![
        permeability_x_1_result,
        permeability_x_2_result,
        permeability_x_3_result,
    ]
}

pub fn compute_mass_conservation_x(lattice: &Lattice) -> Vec<PostResult> {
    let mut mass_conservation_results: Vec<PostResult> = Vec::new();
    for i in 0..lattice.nx {
        let mass_conservation_i = lattice
            .nodes
            .par_iter()
            .filter(|node| node.index[0] == i && node.node_type == NodeType::Fluid)
            .map(|node| {
                let rho = node.density;
                let [ux, uy, uz] = node.velocity;
                let u = (ux * ux + uy * uy + uz * uz).sqrt();
                rho * u
            })
            .sum::<Float>();
        let mass_conservation_i_result = PostResult::new(
            format!("slice_{i:03}"),
            format!("slice ({i:03})").to_string(),
            mass_conservation_i,
            None,
        );
        mass_conservation_results.push(mass_conservation_i_result);
    }
    mass_conservation_results
}

pub fn compute_tpms_quantities(lattice: &Lattice) -> Vec<PostResult> {
    let number_of_solid_nodes = lattice.solid_nodes.as_ref().unwrap().len() as Float;
    let number_of_fluid_nodes = lattice.fluid_nodes.as_ref().unwrap().len() as Float;
    let rho_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.density)
        .sum::<Float>();
    let ux_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.velocity[0])
        .sum::<Float>();
    let uy_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.velocity[1])
        .sum::<Float>();
    let uz_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| node.velocity[2])
        .sum::<Float>();
    let u_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| {
            (node.velocity[0] * node.velocity[0]
                + node.velocity[1] * node.velocity[1]
                + node.velocity[2] * node.velocity[2])
                .sqrt()
        })
        .sum::<Float>();
    let u_max = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| {
            (node.velocity[0] * node.velocity[0]
                + node.velocity[1] * node.velocity[1]
                + node.velocity[2] * node.velocity[2])
                .sqrt()
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let case_parameters = lattice.case_parameters.as_ref().unwrap();
    let porosity = number_of_fluid_nodes / (number_of_fluid_nodes + number_of_solid_nodes);
    let rho_mean = rho_sum / number_of_fluid_nodes;
    let ux_mean = ux_sum / number_of_fluid_nodes;
    let uy_mean = uy_sum / number_of_fluid_nodes;
    let uz_mean = uz_sum / number_of_fluid_nodes;
    let u_mean = u_sum / number_of_fluid_nodes;
    let physical_u_mean = u_mean * case_parameters.velocity_conversion_factor;
    let tortuosity_x = u_sum / ux_sum;
    let tortuosity_y = u_sum / uy_sum;
    let tortuosity_z = u_sum / uz_sum;
    let pressure_inlet_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .filter(|node| node.index[0] == 0)
        .map(|node| {
            let density_prime = node.density - LATTICE_DENSITY;
            let pressure_prime = CS_2 * density_prime;
            let physical_pressure = case_parameters.reference_pressure
                + pressure_prime * case_parameters.pressure_conversion_factor;
            physical_pressure
        })
        .sum::<Float>();
    let pressure_outlet_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .filter(|node| node.index[0] == lattice.nx - 1)
        .map(|node| {
            let density_prime = node.density - LATTICE_DENSITY;
            let pressure_prime = CS_2 * density_prime;
            let physical_pressure = case_parameters.reference_pressure
                + pressure_prime * case_parameters.pressure_conversion_factor;
            physical_pressure
        })
        .sum::<Float>();
    let number_of_liquid_nodes_inlet = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .filter(|node| node.index[0] == 0)
        .count() as Float;
    let number_of_liquid_nodes_outlet = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .filter(|node| node.index[0] == lattice.nx - 1)
        .count() as Float;
    let pressure_inlet_mean = pressure_inlet_sum / number_of_liquid_nodes_inlet;
    let pressure_outlet_mean = pressure_outlet_sum / number_of_liquid_nodes_outlet;
    let pressure_diff_mean = pressure_inlet_mean - pressure_outlet_mean;
    let physical_dynamic_viscosity =
        case_parameters.physical_viscosity * case_parameters.physical_density;
    let permeability_x =
        physical_dynamic_viscosity * lattice.lx * physical_u_mean / pressure_diff_mean;
    let mut pressure_diff_old = 0.0;
    for k in 0..lattice.nz {
        for j in 0..lattice.ny {
            let rho_in = lattice.get_node(&[0, j, k]).density;
            let rho_out = lattice.get_node(&[lattice.nx - 1, j, k]).density;
            let pressure_in =
                CS_2 * (rho_in - LATTICE_DENSITY) + case_parameters.reference_pressure;
            let pressure_out =
                CS_2 * (rho_out - LATTICE_DENSITY) + case_parameters.reference_pressure;
            pressure_diff_old += pressure_in - pressure_out;
        }
    }
    let permeability_x_old =
        physical_dynamic_viscosity * lattice.lx * physical_u_mean / pressure_diff_old;
    let rho_result: PostResult = PostResult::new(
        "mean_density".to_string(),
        "mean density".to_string(),
        rho_mean,
        None,
    );
    let number_of_solid_nodes_result: PostResult = PostResult::new(
        "n_solid_nodes".to_string(),
        "number of solid nodes".to_string(),
        number_of_solid_nodes,
        None,
    );
    let number_of_fluid_nodes_result: PostResult = PostResult::new(
        "n_fluid_nodes".to_string(),
        "number of fluid nodes".to_string(),
        number_of_fluid_nodes,
        None,
    );
    let porosity_result: PostResult = PostResult::new(
        "porosity".to_string(),
        "porosity".to_string(),
        porosity,
        None,
    );
    let u_result: PostResult = PostResult::new(
        "mean_velocity".to_string(),
        "mean velocity (magnitude)".to_string(),
        u_mean,
        None,
    );
    let ux_result: PostResult = PostResult::new(
        "mean_velocity_x".to_string(),
        "mean velocity (x)".to_string(),
        ux_mean,
        None,
    );
    let uy_result: PostResult = PostResult::new(
        "mean_velocity_y".to_string(),
        "mean velocity (y)".to_string(),
        uy_mean,
        None,
    );
    let uz_result: PostResult = PostResult::new(
        "mean_velocity_z".to_string(),
        "mean velocity (z)".to_string(),
        uz_mean,
        None,
    );
    let maximum_velocity_result: PostResult = PostResult::new(
        "max_velocity".to_string(),
        "maximum velocity (magnitude)".to_string(),
        u_max,
        None,
    );
    let tortuosity_x_result: PostResult = PostResult::new(
        "tortuosity_x".to_string(),
        "tortuosity (x)".to_string(),
        tortuosity_x,
        None,
    );
    let tortuosity_y_result: PostResult = PostResult::new(
        "tortuosity_y".to_string(),
        "tortuosity (y)".to_string(),
        tortuosity_y,
        None,
    );
    let tortuosity_z_result: PostResult = PostResult::new(
        "tortuosity_z".to_string(),
        "tortuosity (z)".to_string(),
        tortuosity_z,
        None,
    );
    let pressure_inlet_mean_result: PostResult = PostResult::new(
        "pressure_inlet".to_string(),
        "pressure inlet mean".to_string(),
        pressure_inlet_mean,
        None,
    );
    let pressure_outlet_mean_result: PostResult = PostResult::new(
        "pressure_outlet".to_string(),
        "pressure outlet mean".to_string(),
        pressure_outlet_mean,
        None,
    );
    let pressure_diff_mean_result: PostResult = PostResult::new(
        "pressure_diff".to_string(),
        "pressure difference mean".to_string(),
        pressure_diff_mean,
        None,
    );
    let permeability_x_result: PostResult = PostResult::new(
        "permeability_x".to_string(),
        "permeability (x)".to_string(),
        permeability_x,
        None,
    );
    let pressure_diff_old_result: PostResult = PostResult::new(
        "pres_diff_old".to_string(),
        "pressure difference old".to_string(),
        pressure_diff_old,
        None,
    );
    let permeability_x_old_result: PostResult = PostResult::new(
        "perm_x_old".to_string(),
        "permeability (x) old".to_string(),
        permeability_x_old,
        None,
    );
    vec![
        number_of_solid_nodes_result,
        number_of_fluid_nodes_result,
        porosity_result,
        rho_result,
        u_result,
        ux_result,
        uy_result,
        uz_result,
        maximum_velocity_result,
        tortuosity_x_result,
        tortuosity_y_result,
        tortuosity_z_result,
        pressure_inlet_mean_result,
        pressure_outlet_mean_result,
        pressure_diff_mean_result,
        permeability_x_result,
        pressure_diff_old_result,
        permeability_x_old_result,
    ]
}
