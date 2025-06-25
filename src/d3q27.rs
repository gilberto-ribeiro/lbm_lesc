pub mod bc;
pub mod io;
pub mod post;

pub use bc::{BoundaryCondition, BoundaryFace};

use crate::global_variables::*;
use crate::io::WriteDataMode;
use crate::{NodeType, Residuals};
use rayon::prelude::*;
use std::collections::HashMap;
use std::process;
use std::time::Instant;

pub const D: usize = 3;

pub const Q: usize = 27;

pub const C: [[i32; D]; Q] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, -1],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, -1],
    [0, -1, 1],
    [1, 1, 1],
    [-1, -1, -1],
    [1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [-1, 1, -1],
    [-1, 1, 1],
    [1, -1, -1],
];

pub const W: [Float; Q] = [
    8.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
];

#[derive(Clone)]
pub struct Simulation {
    pub case_name: String,
    pub time_step: usize,
    pub simulation_time: Instant,
    pub tolerance_density: Float,
    pub tolerance_velocity_x: Float,
    pub tolerance_velocity_y: Float,
    pub tolerance_velocity_z: Float,
    pub min_iter: usize,
    pub max_iter: usize,
    pub write_data_mode: WriteDataMode,
}

impl Simulation {
    pub fn next_step(&mut self) {
        self.time_step += 1;
    }

    pub fn stop_condition(&self, residuals: &Residuals) -> bool {
        let converged_density = residuals.density <= self.tolerance_density;
        let converged_velocity_x = residuals.velocity[0] <= self.tolerance_velocity_x;
        let converged_velocity_y = residuals.velocity[1] <= self.tolerance_velocity_y;
        let converged_velocity_z = residuals.velocity[2] <= self.tolerance_velocity_z;
        let converged_quantities = converged_density
            && converged_velocity_x
            && converged_velocity_y
            && converged_velocity_z;
        let min_iterations = self.time_step > self.min_iter;
        let max_iterations = self.time_step > self.max_iter;
        (min_iterations && converged_quantities) || max_iterations
    }
}

impl Simulation {
    pub fn new() -> Self {
        Self {
            case_name: String::from(CASE_NAME),
            time_step: 0,
            simulation_time: Instant::now(),
            tolerance_density: TOLERANCE_DENSITY,
            tolerance_velocity_x: TOLERANCE_VELOCITY_X,
            tolerance_velocity_y: TOLERANCE_VELOCITY_Y,
            tolerance_velocity_z: TOLERANCE_VELOCITY_Z,
            min_iter: MIN_ITER,
            max_iter: MAX_ITER,
            write_data_mode: WriteDataMode::Frequency(100),
        }
    }

    pub fn from_setup(parameters: HashMap<String, String>) -> Self {
        let case_name = parameters["case_name"].clone();
        let tolerance_density = parameters["tolerance_density"].parse::<Float>().unwrap();
        let tolerance_velocity_x = parameters["tolerance_velocity_x"].parse::<Float>().unwrap();
        let tolerance_velocity_y = parameters["tolerance_velocity_y"].parse::<Float>().unwrap();
        let tolerance_velocity_z = parameters["tolerance_velocity_z"].parse::<Float>().unwrap();
        let min_iter = parameters["min_iter"].parse::<usize>().unwrap();
        let max_iter = parameters["max_iter"].parse::<usize>().unwrap();
        let write_data_mode = Simulation::set_write_data_mode(&parameters["write_data_mode"]);
        Self {
            case_name,
            time_step: 0,
            simulation_time: Instant::now(),
            tolerance_density,
            tolerance_velocity_x,
            tolerance_velocity_y,
            tolerance_velocity_z,
            min_iter,
            max_iter,
            write_data_mode,
        }
    }

    fn set_write_data_mode(mode: &str) -> WriteDataMode {
        let mut write_data_mode = mode.split_whitespace();
        let mode = write_data_mode.next().unwrap();
        match mode {
            "frequency" => {
                let frequency = write_data_mode.next().unwrap().parse::<usize>().unwrap();
                WriteDataMode::Frequency(frequency)
            }
            "list" => {
                let list = write_data_mode
                    .map(|x| x.parse::<usize>().unwrap())
                    .collect();
                WriteDataMode::ListOfSteps(list)
            }
            _ => panic!("Invalid write results mode: {mode}."),
        }
    }
}

#[derive(Clone)]
pub struct Node {
    pub node_type: NodeType,

    pub boundary_faces: Option<Vec<BoundaryFace>>,

    pub index: [usize; D],

    pub coordinates: [Float; D],

    pub density: Float,

    pub velocity: [Float; D],

    pub f: Box<[Float; Q]>,

    pub f_eq: Box<[Float; Q]>,

    pub f_star: Box<[Float; Q]>,
}

impl Node {
    pub fn update_density(&mut self) {
        match (IGNORE_SOLID, self.node_type) {
            (true, NodeType::Solid) => {}
            _ => {
                self.density = self.f.iter().sum();
            }
        }
    }

    pub fn update_velocity(&mut self) {
        match (IGNORE_SOLID, self.node_type) {
            (true, NodeType::Solid) => {}
            _ => {
                self.velocity[0] = (1.0 / self.density)
                    * (self.f[1] - self.f[2] + self.f[7] - self.f[8] + self.f[9] - self.f[10]
                        + self.f[13]
                        - self.f[14]
                        + self.f[15]
                        - self.f[16]
                        + self.f[19]
                        - self.f[20]
                        + self.f[21]
                        - self.f[22]
                        + self.f[23]
                        - self.f[24]
                        - self.f[25]
                        + self.f[26]);
                self.velocity[1] = (1.0 / self.density)
                    * (self.f[3] - self.f[4] + self.f[7] - self.f[8] + self.f[11]
                        - self.f[12]
                        - self.f[13]
                        + self.f[14]
                        + self.f[17]
                        - self.f[18]
                        + self.f[19]
                        - self.f[20]
                        + self.f[21]
                        - self.f[22]
                        - self.f[23]
                        + self.f[24]
                        + self.f[25]
                        - self.f[26]);
                self.velocity[2] = (1.0 / self.density)
                    * (self.f[5] - self.f[6] + self.f[9] - self.f[10] + self.f[11]
                        - self.f[12]
                        - self.f[15]
                        + self.f[16]
                        - self.f[17]
                        + self.f[18]
                        + self.f[19]
                        - self.f[20]
                        - self.f[21]
                        + self.f[22]
                        + self.f[23]
                        - self.f[24]
                        + self.f[25]
                        - self.f[26]);
            }
        }
    }

    pub fn update_moments(&mut self) {
        self.update_density();
        self.update_velocity();
    }

    pub fn equilibrium(&mut self, explicit_version: bool) {
        match (IGNORE_SOLID, self.node_type) {
            (true, NodeType::Solid) => {}
            _ => {
                let [ux, uy, uz] = self.velocity;
                let ux_2 = ux * ux;
                let uy_2 = uy * uy;
                let uz_2 = uz * uz;
                let u_2 = ux_2 + uy_2 + uz_2;
                if explicit_version {
                    todo!();
                } else {
                    for q in 0..Q {
                        let cx = C[q][0] as Float;
                        let cy = C[q][1] as Float;
                        let cz = C[q][2] as Float;
                        let u_dot_c = ux * cx + uy * cy + uz * cz;
                        self.f_eq[q] = W[q]
                            * self.density
                            * (1.0 + CS_2_INV * u_dot_c + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                - 0.5 * CS_2_INV * u_2);
                    }
                }
            }
        }
    }

    pub fn collision_step(&mut self, omega: Float, omega_prime: Float) {
        match (IGNORE_SOLID, self.node_type) {
            (true, NodeType::Solid) => {}
            _ => {
                for q in 0..Q {
                    self.f_star[q] = omega_prime * self.f[q] + omega * self.f_eq[q];
                }
            }
        }
    }
}

impl Node {
    pub fn new(density: Float, velocity: [Float; D], omega: Float, omega_prime: Float) -> Self {
        let mut node = Self {
            node_type: NodeType::Fluid,
            boundary_faces: None,
            index: [0; D],
            coordinates: [0.0; D],
            density,
            velocity,
            f: Box::new([0.0; Q]),
            f_eq: Box::new([0.0; Q]),
            f_star: Box::new([0.0; Q]),
        };
        node.equilibrium(false);
        node.f = node.f_eq.clone();
        node.collision_step(omega, omega_prime);
        node
    }
}

#[derive(Copy, Clone)]
pub struct ShallowNode {
    pub coordinates: [Float; D],

    pub density: Float,

    pub velocity: [Float; D],
}

impl ShallowNode {
    pub fn new(density: Float, velocity: [Float; D]) -> Self {
        Self {
            coordinates: [0.0; D],
            density,
            velocity,
        }
    }
}

pub struct Lattice {
    pub dimensions: usize,
    pub delta_x: Float,
    pub delta_t: Float,
    pub tau: Float,
    pub omega: Float,
    pub omega_prime: Float,
    pub physical_density: Float,
    pub reference_pressure: Float,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub lx: Float,
    pub ly: Float,
    pub lz: Float,
    pub nodes: Vec<Node>,
    pub bounce_back_map: bool,
    pub boundary_conditions: HashMap<BoundaryFace, BoundaryCondition>,
    pub solid_nodes: Option<Vec<[usize; D]>>,
    pub fluid_nodes: Option<Vec<[usize; D]>>,
    pub interior_nodes: Option<Vec<[usize; D]>>,
    pub case_parameters: Option<post::vtk::CaseParameters>,
}

impl Lattice {
    pub fn get_node(&self, index: &[usize]) -> &Node {
        let i = index[0];
        let j = index[1];
        let k = index[2];
        let number_of_slices = self.nx * self.ny * k;
        let number_of_rows = self.nx * j;
        let index = i + number_of_rows + number_of_slices;
        &self.nodes[index]
    }

    pub fn get_node_mut(&mut self, index: &[usize]) -> &mut Node {
        let i = index[0];
        let j = index[1];
        let k = index[2];
        let number_of_slices = self.nx * self.ny * k;
        let number_of_rows = self.nx * j;
        let index = i + number_of_rows + number_of_slices;
        &mut self.nodes[index]
    }

    pub fn update_moments(&mut self) {
        self.nodes
            .par_iter_mut()
            .for_each(|node| node.update_moments());
    }

    pub fn equilibrium(&mut self, explicit_version: bool) {
        self.nodes
            .par_iter_mut()
            .for_each(|node| node.equilibrium(explicit_version));
    }

    pub fn collision_step(&mut self) {
        self.nodes
            .par_iter_mut()
            .for_each(|node| node.collision_step(self.omega, self.omega_prime));
    }

    pub fn streaming_step(&mut self) {
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
                    for q in 0..Q {
                        let [cx, cy, cz] = C[q];
                        let new_i = ((i as i32) + cx).rem_euclid(self.nx as i32) as usize;
                        let new_j = ((j as i32) + cy).rem_euclid(self.ny as i32) as usize;
                        let new_k = ((k as i32) + cz).rem_euclid(self.nz as i32) as usize;
                        self.get_node_mut(&[new_i, new_j, new_k]).f[q] =
                            self.get_node(&[i, j, k]).f_star[q];
                    }
                }
            }
        }
    }

    pub fn boundary_condition(&mut self) {
        let boundary_conditions: Vec<_> = self
            .boundary_conditions
            .iter()
            .map(|(&k, v)| (k, v.clone()))
            .collect();
        for (boundary_face, boundary_condition) in boundary_conditions {
            self.apply_boundary_condition(boundary_face, boundary_condition);
        }
    }

    pub fn initialize_bounce_back_map(&mut self, code: i32) {
        if self.bounce_back_map {
            let map = io::read_bounce_back_map().unwrap();
            self.apply_bounce_back_map(map, code);
        }
    }

    pub fn apply_bounce_back_map(&mut self, map: Vec<Vec<Vec<i32>>>, code: i32) {
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
                    if map[k][j][i] == code {
                        self.get_node_mut(&[i, j, k]).node_type = NodeType::Solid;
                    }
                }
            }
        }
        self.set_node_type();
    }

    pub fn par_bounce_back(&mut self, bounce_back_matrix: Vec<Vec<Vec<NodeType>>>) {
        if self.bounce_back_map {
            self.nodes
                .par_iter_mut()
                .filter(|node| self.interior_nodes.as_ref().unwrap().contains(&node.index))
                .for_each(|node| {
                    for q in 0..Q {
                        let [cx, cy, cz] = C[q];
                        let [i, j, k] = node.index;
                        let new_i = ((i as i32) + cx).rem_euclid(self.nx as i32) as usize;
                        let new_j = ((j as i32) + cy).rem_euclid(self.ny as i32) as usize;
                        let new_k = ((k as i32) + cz).rem_euclid(self.nz as i32) as usize;
                        if bounce_back_matrix[new_k][new_j][new_i] == NodeType::Solid {
                            node.f[q] = node.f_star[q];
                        }
                    }
                });
        }
    }

    pub fn bounce_back(&mut self) {
        if self.bounce_back_map {
            for i in 1..self.nx - 1 {
                for j in 1..self.ny - 1 {
                    for k in 1..self.nz - 1 {
                        for q in 0..Q {
                            let [cx, cy, cz] = C[q];
                            let new_i = ((i as i32) + cx) as usize;
                            let new_j = ((j as i32) + cy) as usize;
                            let new_k = ((k as i32) + cz) as usize;
                            if self.get_node(&[new_i, new_j, new_k]).node_type == NodeType::Solid {
                                self.get_node_mut(&[i, j, k]).f[q] =
                                    self.get_node(&[i, j, k]).f_star[q];
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn generate_shallow_lattice(&self) -> ShallowLattice {
        let mut shallow_lattice = ShallowLattice::new(self.nx, self.ny, self.nz, 1.0, [0.0; D]);
        shallow_lattice.update_shallow_lattice(&self);
        shallow_lattice
    }

    pub fn generate_bounce_back_matrix(&self) -> Vec<Vec<Vec<NodeType>>> {
        let mut bounce_back_matrix = vec![vec![vec![NodeType::Fluid; self.nx]; self.ny]; self.nz];
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    bounce_back_matrix[k][j][i] = self.get_node(&[i, j, k]).node_type;
                }
            }
        }
        bounce_back_matrix
    }

    pub fn compute_residuals(&self, old_lattice: &ShallowLattice) -> Residuals {
        let density = self
            .nodes
            .par_iter()
            .zip(old_lattice.nodes.par_iter())
            .map(|(node, old_node)| (node.density - old_node.density).powi(2))
            .sum::<Float>()
            .sqrt();
        let velocity_x = self
            .nodes
            .par_iter()
            .zip(old_lattice.nodes.par_iter())
            .map(|(node, old_node)| (node.velocity[0] - old_node.velocity[0]).powi(2))
            .sum::<Float>()
            .sqrt();
        let velocity_y = self
            .nodes
            .par_iter()
            .zip(old_lattice.nodes.par_iter())
            .map(|(node, old_node)| (node.velocity[1] - old_node.velocity[1]).powi(2))
            .sum::<Float>()
            .sqrt();
        let velocity_z = self
            .nodes
            .par_iter()
            .zip(old_lattice.nodes.par_iter())
            .map(|(node, old_node)| (node.velocity[2] - old_node.velocity[2]).powi(2))
            .sum::<Float>()
            .sqrt();
        Residuals {
            density,
            velocity: vec![velocity_x, velocity_y, velocity_z],
        }
    }

    fn set_node_type(&mut self) {
        let mut solid_nodes = Vec::new();
        let mut fluid_nodes = Vec::new();
        self.nodes.iter().for_each(|node| match node.node_type {
            NodeType::Solid => solid_nodes.push(node.index),
            NodeType::Fluid => fluid_nodes.push(node.index),
        });
        self.solid_nodes = Some(solid_nodes);
        self.fluid_nodes = Some(fluid_nodes);
    }

    fn set_node_boundary_faces(&mut self) {
        self.nodes.par_iter_mut().for_each(|node| {
            let mut boundary_faces = Vec::new();
            let [i, j, k] = node.index;
            if i == 0 {
                boundary_faces.push(BoundaryFace::West);
            }
            if i == self.nx - 1 {
                boundary_faces.push(BoundaryFace::East);
            }
            if j == 0 {
                boundary_faces.push(BoundaryFace::South);
            }
            if j == self.ny - 1 {
                boundary_faces.push(BoundaryFace::North);
            }
            if k == 0 {
                boundary_faces.push(BoundaryFace::Bottom);
            }
            if k == self.nz - 1 {
                boundary_faces.push(BoundaryFace::Top);
            }
            node.boundary_faces = Some(boundary_faces);
        });
    }

    fn set_node_index(&mut self) {
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    self.get_node_mut(&[i, j, k]).index = [i, j, k];
                }
            }
        }
    }

    fn set_node_coordinates(&mut self) {
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    self.get_node_mut(&[i, j, k]).coordinates = [
                        (i as Float) * self.delta_x + 0.5 * self.delta_x,
                        (j as Float) * self.delta_x + 0.5 * self.delta_x,
                        (k as Float) * self.delta_x + 0.5 * self.delta_x,
                    ];
                }
            }
        }
    }

    fn set_interior_nodes(&mut self) {
        let mut interior_nodes = Vec::new();
        for k in 1..self.nz - 1 {
            for j in 1..self.ny - 1 {
                for i in 1..self.nx - 1 {
                    interior_nodes.push([i, j, k]);
                }
            }
        }
        self.interior_nodes = Some(interior_nodes);
    }
}

impl Lattice {
    pub fn initialization(conditions: HashMap<String, String>) -> Self {
        let delta_x = conditions["delta_x"].parse::<Float>().unwrap();
        let delta_t = conditions["delta_t"].parse::<Float>().unwrap();
        let tau = conditions["tau"].parse::<Float>().unwrap();
        let omega = DELTA_T / tau;
        let omega_prime = 1.0 - omega;
        let physical_density = conditions["physical_density"].parse::<Float>().unwrap();
        let reference_pressure = conditions["reference_pressure"].parse::<Float>().unwrap();
        let lx = conditions["lx"].parse::<Float>().unwrap();
        let ly = conditions["ly"].parse::<Float>().unwrap();
        let lz = conditions["lz"].parse::<Float>().unwrap();
        let nx = (lx / delta_x) as usize;
        let ny = (ly / delta_x) as usize;
        let nz = (lz / delta_x) as usize;
        let density = conditions["initial_density"].parse::<Float>().unwrap();
        let velocity = conditions["initial_velocity"]
            .split_whitespace()
            .map(|x| x.parse::<Float>().unwrap())
            .collect::<Vec<Float>>()
            .try_into()
            .unwrap();
        let bounce_back_map = conditions["bounce_back_map"].parse::<bool>().unwrap();
        let boundary_conditions = Self::set_boundary_conditions(conditions);
        let mut lattice = Self {
            dimensions: D,
            delta_x,
            delta_t,
            tau,
            omega,
            omega_prime,
            physical_density,
            reference_pressure,
            nx,
            ny,
            nz,
            lx,
            ly,
            lz,
            nodes: vec![Node::new(density, velocity, omega, omega_prime); nx * ny * nz],
            bounce_back_map,
            boundary_conditions,
            solid_nodes: None,
            fluid_nodes: None,
            interior_nodes: None,
            case_parameters: None,
        };
        let case_parameters = post::vtk::CaseParameters::new(&lattice);
        lattice.case_parameters = Some(case_parameters);
        lattice.set_node_index();
        lattice.set_node_coordinates();
        lattice.set_node_type();
        lattice.set_node_boundary_faces();
        lattice.set_interior_nodes();
        lattice
    }

    fn set_boundary_conditions(
        conditions: HashMap<String, String>,
    ) -> HashMap<BoundaryFace, BoundaryCondition> {
        let conditions = conditions
            .iter()
            .filter(|(k, _)| k.ends_with("_boundary_condition"));
        let mut boundary_conditions = HashMap::new();
        for (k, v) in conditions {
            let face = k.split("_").next().unwrap();
            let boundary_face = match face {
                "north" => BoundaryFace::North,
                "south" => BoundaryFace::South,
                "east" => BoundaryFace::East,
                "west" => BoundaryFace::West,
                "top" => BoundaryFace::Top,
                "bottom" => BoundaryFace::Bottom,
                _ => panic!("Invalid boundary face: {face}."),
            };
            let condition = v.split_whitespace().collect::<Vec<&str>>();
            let boundary_condition = match condition[0] {
                "no_slip" => BoundaryCondition::NoSlip,
                "periodic" => BoundaryCondition::Periodic,
                "dirichlet" => BoundaryCondition::Dirichlet {
                    density: condition[1].parse::<Float>().unwrap(),
                    velocity: [
                        condition[2].parse::<Float>().unwrap(),
                        condition[3].parse::<Float>().unwrap(),
                        condition[4].parse::<Float>().unwrap(),
                    ],
                },
                "pressure_outlet" => BoundaryCondition::PressureOutlet {
                    density: condition[1].parse::<Float>().unwrap(),
                },
                _ => panic!("Invalid boundary condition."),
            };
            boundary_conditions.insert(boundary_face, boundary_condition);
        }
        boundary_conditions
    }
}

pub struct ShallowLattice {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub nodes: Vec<ShallowNode>,
}

impl ShallowLattice {
    pub fn get_node(&self, index: &[usize]) -> &ShallowNode {
        let i = index[0];
        let j = index[1];
        let k = index[2];
        let number_of_slices = self.nx * self.ny * k;
        let number_of_rows = self.nx * j;
        let index = i + number_of_rows + number_of_slices;
        &self.nodes[index]
    }

    pub fn get_node_mut(&mut self, index: &[usize]) -> &mut ShallowNode {
        let i = index[0];
        let j = index[1];
        let k = index[2];
        let number_of_slices = self.nx * self.ny * k;
        let number_of_rows = self.nx * j;
        let index = i + number_of_rows + number_of_slices;
        &mut self.nodes[index]
    }

    pub fn update_shallow_lattice(&mut self, lattice: &Lattice) {
        self.nodes
            .par_iter_mut()
            .zip(lattice.nodes.par_iter())
            .for_each(|(shallow_node, node)| {
                shallow_node.coordinates = node.coordinates;
                shallow_node.density = node.density;
                shallow_node.velocity = node.velocity;
            });
    }
}

impl ShallowLattice {
    pub fn new(nx: usize, ny: usize, nz: usize, density: Float, velocity: [Float; D]) -> Self {
        Self {
            nx,
            ny,
            nz,
            nodes: vec![ShallowNode::new(density, velocity); nx * ny * nz],
        }
    }
}

pub fn run() {
    let mut handle_write_data;

    let mut simulation = Simulation::build_case_setup().unwrap();

    let mut lattice = Lattice::build_case_conditions().unwrap();

    lattice.initialize_bounce_back_map(1);

    if let Err(e) = simulation.write_node_type_vtk(&lattice) {
        eprintln!("Error while writing node type vtk file: {e}.");
        process::exit(1);
    }

    if let Err(e) = simulation.write_post_processing_from_each_n_steps(
        &lattice,
        1,
        post::compute_porosity,
        "porosity.dat",
    ) {
        eprintln!("Error while writing the porosity file: {e}.");
        process::exit(1);
    };

    loop {
        let mut old_lattice = lattice.generate_shallow_lattice();

        lattice.update_moments();
        lattice.equilibrium(false);
        lattice.collision_step();

        lattice.streaming_step();

        lattice.boundary_condition();

        lattice.bounce_back();

        let residuals = lattice.compute_residuals(&old_lattice);

        old_lattice.update_shallow_lattice(&lattice);
        let simulation_clone = simulation.clone();
        handle_write_data = std::thread::spawn(move || simulation_clone.write_data(&old_lattice));

        if let Err(e) = simulation.write_post_processing_from_each_n_steps(
            &lattice,
            1,
            post::compute_mean_density,
            "mean_density.dat",
        ) {
            eprintln!("Error while writing the mean_density file: {e}.");
            process::exit(1);
        };

        if let Err(e) = simulation.write_post_processing_from_each_n_steps(
            &lattice,
            1,
            post::compute_mean_velocities,
            "mean_velocities.dat",
        ) {
            eprintln!("Error while writing the mean_velocities file: {e}.");
            process::exit(1);
        };

        if let Err(e) = simulation.write_post_processing_from_each_n_steps(
            &lattice,
            1,
            post::compute_tortuosities,
            "tortuosities.dat",
        ) {
            eprintln!("Error while writing the tortuosities file: {e}.");
            process::exit(1);
        };

        if let Err(e) = simulation.write_post_processing_from_each_n_steps(
            &lattice,
            1,
            post::compute_max_velocity,
            "max_velocity.dat",
        ) {
            eprintln!("Error while writing the max velocity file: {e}.");
            process::exit(1);
        };

        if let Err(e) = simulation.write_post_processing_from_each_n_steps(
            &lattice,
            1,
            post::compute_mean_pressures_x,
            "mean_pressures_x.dat",
        ) {
            eprintln!("Error while writing the mean pressures x file: {e}.");
            process::exit(1);
        };

        if simulation.stop_condition(&residuals) {
            let old_lattice = lattice.generate_shallow_lattice();
            simulation.write_data_from_steps(&old_lattice);
            simulation.write_vtk_from_steps(&old_lattice);
            break;
        }

        simulation.print_residuals(&residuals);
        if let Err(e) = simulation.write_residuals(&residuals) {
            eprintln!("Error while writing the residuals file: {e}.");
            process::exit(1);
        };

        simulation.next_step();
    }

    handle_write_data.join().unwrap();
}

pub fn run_benchmark() {
    let mut handle_write_results;

    let bcs_time = Instant::now();
    let mut simulation = Simulation::build_case_setup().unwrap();
    let bcs_duration = bcs_time.elapsed();

    let bcc_time = Instant::now();
    let mut lattice = Lattice::build_case_conditions().unwrap();
    let bcc_duration = bcc_time.elapsed();

    let ibbm_time = Instant::now();
    lattice.initialize_bounce_back_map(1);
    let ibbm_duration = ibbm_time.elapsed();

    let wnt_time = Instant::now();
    if let Err(e) = simulation.write_node_type_vtk(&lattice) {
        eprintln!("Error while writing node type vtk file: {e}.");
        process::exit(1);
    }
    let wnt_duration = wnt_time.elapsed();

    if let Err(e) = simulation.write_post_processing_from_each_n_steps(
        &lattice,
        1,
        post::compute_porosity,
        "porosity.dat",
    ) {
        eprintln!("Error while writing the porosity file: {e}.");
        process::exit(1);
    };

    loop {
        let loop_time = Instant::now();

        let gsl_time = Instant::now();
        let mut old_lattice = lattice.generate_shallow_lattice();
        let gsl_duration = gsl_time.elapsed();

        let um_time = Instant::now();
        lattice.update_moments();
        let um_duration = um_time.elapsed();

        let e_time = Instant::now();
        lattice.equilibrium(false);
        let e_duration = e_time.elapsed();

        let cs_time = Instant::now();
        lattice.collision_step();
        let cs_duration = cs_time.elapsed();

        let ss_time = Instant::now();
        lattice.streaming_step();
        let ss_duration = ss_time.elapsed();

        let bc_time = Instant::now();
        lattice.boundary_condition();
        let bc_duration = bc_time.elapsed();

        let bb_time = Instant::now();

        lattice.bounce_back();
        let bb_duration = bb_time.elapsed();

        let cr_time = Instant::now();
        let residuals = lattice.compute_residuals(&old_lattice);
        let cr_duration = cr_time.elapsed();

        let wr_time = Instant::now();
        old_lattice.update_shallow_lattice(&lattice);
        let simulation_clone = simulation.clone();
        handle_write_results =
            std::thread::spawn(move || simulation_clone.write_data(&old_lattice));
        let wr_duration = wr_time.elapsed();

        let cmd_time = Instant::now();
        if let Err(e) = simulation.write_post_processing_from_each_n_steps(
            &lattice,
            1,
            post::compute_mean_density,
            "mean_density.dat",
        ) {
            eprintln!("Error while writing the mean_density file: {e}.");
            process::exit(1);
        };
        let cmd_duration = cmd_time.elapsed();

        let cmv_time = Instant::now();
        if let Err(e) = simulation.write_post_processing_from_each_n_steps(
            &lattice,
            1,
            post::compute_mean_velocities,
            "mean_velocities.dat",
        ) {
            eprintln!("Error while writing the mean_velocities file: {e}.");
            process::exit(1);
        };
        let cmv_duration = cmv_time.elapsed();

        if let Err(e) = simulation.write_post_processing_from_each_n_steps(
            &lattice,
            1,
            post::compute_tortuosities,
            "tortuosities.dat",
        ) {
            eprintln!("Error while writing the tortuosities file: {e}.");
            process::exit(1);
        };

        let sc_time = Instant::now();
        if simulation.stop_condition(&residuals) {
            let old_lattice = lattice.generate_shallow_lattice();
            simulation.write_data_from_steps(&old_lattice);
            simulation.write_vtk_from_steps(&old_lattice);
            break;
        }
        let sc_duration = sc_time.elapsed();

        let pr_time = Instant::now();
        simulation.print_residuals(&residuals);
        if let Err(e) = simulation.write_residuals(&residuals) {
            eprintln!("Error while writing the residuals file: {e}.");
            process::exit(1);
        };
        let pr_duration = pr_time.elapsed();

        let loop_duration = loop_time.elapsed();

        let elapsed_times = [
            ("bcs", bcs_duration),
            ("bcc", bcc_duration),
            ("ibbm", ibbm_duration),
            ("wnt", wnt_duration),
            ("gsl", gsl_duration),
            ("um", um_duration),
            ("e", e_duration),
            ("cs", cs_duration),
            ("ss", ss_duration),
            ("bc", bc_duration),
            ("bb", bb_duration),
            ("cr", cr_duration),
            ("wr", wr_duration),
            ("cmd", cmd_duration),
            ("cmv", cmv_duration),
            ("sc", sc_duration),
            ("pr", pr_duration),
            ("loop", loop_duration),
        ];

        let _ = crate::io::write_inside_loop_elapsed_time(&elapsed_times, &simulation.time_step);

        simulation.next_step();
    }

    handle_write_results.join().unwrap();
}
