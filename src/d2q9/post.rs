pub mod vtk;

use super::*;
use crate::global_variables::*;
use crate::post::PostResult;
use rayon::prelude::*;

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
    let u_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| {
            (node.velocity[0] * node.velocity[0] + node.velocity[1] * node.velocity[1]).sqrt()
        })
        .sum::<Float>();
    let number_of_fluid_nodes = lattice.fluid_nodes.as_ref().unwrap().len() as Float;
    let ux_mean = ux_sum / number_of_fluid_nodes;
    let uy_mean = uy_sum / number_of_fluid_nodes;
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
    vec![u_result, ux_result, uy_result]
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
    let u_sum = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| {
            (node.velocity[0] * node.velocity[0] + node.velocity[1] * node.velocity[1]).sqrt()
        })
        .sum::<Float>();
    let tortuosity_x = u_sum / ux_sum;
    let tortuosity_y = u_sum / uy_sum;
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
    vec![tortuosity_x_result, tortuosity_y_result]
}

pub fn compute_max_velocity(lattice: &Lattice) -> Vec<PostResult> {
    let max_velocity = lattice
        .nodes
        .par_iter()
        .filter(|node| node.node_type == NodeType::Fluid)
        .map(|node| {
            (node.velocity[0] * node.velocity[0] + node.velocity[1] * node.velocity[1]).sqrt()
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
