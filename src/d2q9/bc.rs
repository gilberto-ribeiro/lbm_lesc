use super::*;
use crate::global_variables::*;
use rayon::prelude::*;

const Q_BAR: [usize; Q] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

const Q_EAST: [usize; 3] = [1, 5, 8];

const Q_WEST: [usize; 3] = [3, 6, 7];

const Q_NORTH: [usize; 3] = [2, 5, 6];

const Q_SOUTH: [usize; 3] = [4, 7, 8];

const ZOUHE_ERROR_MESSAGE: &'static str =
    "Error: Wrong input of Zou & He (1997) boundary conditions.";

#[derive(Clone)]
pub enum BoundaryCondition {
    NoSlip,
    Dirichlet {
        density: Float,
        velocity: [Float; D],
    },
    PressureOutlet {
        density: Float,
    },
    ZouHe {
        density: Option<Float>,
        velocity: [Option<Float>; D],
    },
    Periodic,
}

#[derive(Hash, Eq, PartialEq, Clone, Copy)]
pub enum BoundaryFace {
    East,
    West,
    North,
    South,
}

impl Lattice {
    pub fn apply_boundary_condition(
        &mut self,
        boundary_face: BoundaryFace,
        boundary_condition: BoundaryCondition,
    ) {
        match boundary_condition {
            BoundaryCondition::NoSlip => self.no_slip(boundary_face),
            BoundaryCondition::Dirichlet { density, velocity } => {
                self.dirichlet(boundary_face, density, velocity)
            }
            BoundaryCondition::PressureOutlet { density } => {
                self.pressure_outlet(boundary_face, density)
            }
            BoundaryCondition::ZouHe { density, velocity } => {
                self.zou_he(boundary_face, density, velocity)
            }
            BoundaryCondition::Periodic => {}
        }
    }

    fn no_slip(&mut self, boundary_face: BoundaryFace) {
        let q_face = match boundary_face {
            BoundaryFace::East => Q_EAST,
            BoundaryFace::West => Q_WEST,
            BoundaryFace::North => Q_NORTH,
            BoundaryFace::South => Q_SOUTH,
        };
        self.nodes
            .par_iter_mut()
            .filter(|node| {
                node.boundary_faces
                    .as_ref()
                    .unwrap()
                    .contains(&boundary_face)
            })
            .for_each(|node| {
                inner_no_slip(node, q_face);
            });
    }

    fn dirichlet(&mut self, boundary_face: BoundaryFace, density: Float, velocity: [Float; D]) {
        let q_face = match boundary_face {
            BoundaryFace::East => Q_EAST,
            BoundaryFace::West => Q_WEST,
            BoundaryFace::North => Q_NORTH,
            BoundaryFace::South => Q_SOUTH,
        };
        self.nodes
            .par_iter_mut()
            .filter(|node| {
                node.boundary_faces
                    .as_ref()
                    .unwrap()
                    .contains(&boundary_face)
            })
            .for_each(|node| match node.node_type {
                NodeType::Fluid => {
                    for q in q_face {
                        let u_dot_c = velocity
                            .iter()
                            .zip(C[q].iter())
                            .map(|(&u, &c)| u * (c as Float))
                            .sum::<Float>();
                        node.f[Q_BAR[q]] =
                            node.f_star[q] - 2.0 * W[q] * density * CS_2_INV * u_dot_c;
                    }
                }
                NodeType::Solid => {
                    inner_no_slip(node, q_face);
                }
            });
    }

    fn pressure_outlet(&mut self, boundary_face: BoundaryFace, density: Float) {
        match boundary_face {
            BoundaryFace::East => {
                for j in 0..self.ny {
                    match self.get_node(&[self.nx - 1, j]).node_type {
                        NodeType::Fluid => {
                            let ux = 1.5 * self.get_node(&[self.nx - 1, j]).velocity[0]
                                - 0.5 * self.get_node(&[self.nx - 2, j]).velocity[0];
                            let uy = 1.5 * self.get_node(&[self.nx - 1, j]).velocity[1]
                                - 0.5 * self.get_node(&[self.nx - 2, j]).velocity[1];
                            let u_2 = ux * ux + uy * uy;
                            for q in Q_EAST {
                                let cx = C[q][0] as Float;
                                let cy = C[q][1] as Float;
                                let u_dot_c = ux * cx + uy * cy;
                                self.get_node_mut(&[self.nx - 1, j]).f[Q_BAR[q]] =
                                    -self.get_node(&[self.nx - 1, j]).f_star[q]
                                        + 2.0
                                            * W[q]
                                            * density
                                            * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                - 0.5 * CS_2_INV * u_2);
                            }
                        }
                        NodeType::Solid => {
                            for q in Q_EAST {
                                self.get_node_mut(&[self.nx - 1, j]).f[Q_BAR[q]] =
                                    self.get_node(&[self.nx - 1, j]).f_star[q];
                            }
                        }
                    }
                }
            }
            BoundaryFace::West => {
                for j in 0..self.ny {
                    match self.get_node(&[0, j]).node_type {
                        NodeType::Fluid => {
                            let ux = 1.5 * self.get_node(&[0, j]).velocity[0]
                                - 0.5 * self.get_node(&[1, j]).velocity[0];
                            let uy = 1.5 * self.get_node(&[0, j]).velocity[1]
                                - 0.5 * self.get_node(&[1, j]).velocity[1];
                            let u_2 = ux * ux + uy * uy;
                            for q in Q_WEST {
                                let cx = C[q][0] as Float;
                                let cy = C[q][1] as Float;
                                let u_dot_c = ux * cx + uy * cy;
                                self.get_node_mut(&[0, j]).f[Q_BAR[q]] =
                                    -self.get_node(&[0, j]).f_star[q]
                                        + 2.0
                                            * W[q]
                                            * density
                                            * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                - 0.5 * CS_2_INV * u_2);
                            }
                        }
                        NodeType::Solid => {
                            for q in Q_WEST {
                                self.get_node_mut(&[0, j]).f[Q_BAR[q]] =
                                    self.get_node(&[0, j]).f_star[q];
                            }
                        }
                    }
                }
            }
            BoundaryFace::North => {
                for i in 0..self.nx {
                    match self.get_node(&[i, self.ny - 1]).node_type {
                        NodeType::Fluid => {
                            let ux = 1.5 * self.get_node(&[i, self.ny - 1]).velocity[0]
                                - 0.5 * self.get_node(&[i, self.ny - 2]).velocity[0];
                            let uy = 1.5 * self.get_node(&[i, self.ny - 1]).velocity[1]
                                - 0.5 * self.get_node(&[i, self.ny - 2]).velocity[1];
                            let u_2 = ux * ux + uy * uy;
                            for q in Q_NORTH {
                                let cx = C[q][0] as Float;
                                let cy = C[q][1] as Float;
                                let u_dot_c = ux * cx + uy * cy;
                                self.get_node_mut(&[i, self.ny - 1]).f[Q_BAR[q]] =
                                    -self.get_node(&[i, self.ny - 1]).f_star[q]
                                        + 2.0
                                            * W[q]
                                            * density
                                            * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                - 0.5 * CS_2_INV * u_2);
                            }
                        }
                        NodeType::Solid => {
                            for q in Q_NORTH {
                                self.get_node_mut(&[i, self.ny - 1]).f[Q_BAR[q]] =
                                    self.get_node(&[i, self.ny - 1]).f_star[q];
                            }
                        }
                    }
                }
            }
            BoundaryFace::South => {
                for i in 0..self.nx {
                    match self.get_node(&[i, 0]).node_type {
                        NodeType::Fluid => {
                            let ux = 1.5 * self.get_node(&[i, 0]).velocity[0]
                                - 0.5 * self.get_node(&[i, 1]).velocity[0];
                            let uy = 1.5 * self.get_node(&[i, 0]).velocity[1]
                                - 0.5 * self.get_node(&[i, 1]).velocity[1];
                            let u_2 = ux * ux + uy * uy;
                            for q in Q_SOUTH {
                                let cx = C[q][0] as Float;
                                let cy = C[q][1] as Float;
                                let u_dot_c = ux * cx + uy * cy;
                                self.get_node_mut(&[i, 0]).f[Q_BAR[q]] =
                                    -self.get_node(&[i, 0]).f_star[q]
                                        + 2.0
                                            * W[q]
                                            * density
                                            * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                - 0.5 * CS_2_INV * u_2);
                            }
                        }
                        NodeType::Solid => {
                            for q in Q_SOUTH {
                                self.get_node_mut(&[i, 0]).f[Q_BAR[q]] =
                                    self.get_node(&[i, 0]).f_star[q];
                            }
                        }
                    }
                }
            }
        }
    }

    fn zou_he(
        &mut self,
        boundary_face: BoundaryFace,
        density: Option<Float>,
        velocity: [Option<Float>; D],
    ) {
        match boundary_face {
            BoundaryFace::East => {
                self.nodes
                    .par_iter_mut()
                    .filter(|node| {
                        node.boundary_faces
                            .as_ref()
                            .unwrap()
                            .contains(&boundary_face)
                    })
                    .for_each(|node| match node.node_type {
                        NodeType::Fluid => match (density, velocity) {
                            (None, [Some(ux), Some(uy)]) => {
                                let f = &mut node.f;
                                let rho = 1.0 / (1.0 + ux)
                                    * (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8]));
                                f[3] = f[1] - (2.0 / 3.0) * rho * ux;
                                f[7] = f[5] + 0.5 * (f[2] - f[4])
                                    - 0.5 * rho * uy
                                    - (1.0 / 6.0) * rho * ux;
                                f[6] = f[8] - 0.5 * (f[2] - f[4]) + 0.5 * rho * uy
                                    - (1.0 / 6.0) * rho * ux;
                            }
                            (Some(rho), [None, Some(uy)]) => {
                                let f = &mut node.f;
                                let ux = 1.0 / rho
                                    * (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8]))
                                    - 1.0;
                                f[3] = f[1] - (2.0 / 3.0) * rho * ux;
                                f[7] = f[5] + 0.5 * (f[2] - f[4])
                                    - 0.5 * rho * uy
                                    - (1.0 / 6.0) * rho * ux;
                                f[6] = f[8] - 0.5 * (f[2] - f[4]) + 0.5 * rho * uy
                                    - (1.0 / 6.0) * rho * ux;
                            }
                            _ => {
                                panic!("{}", ZOUHE_ERROR_MESSAGE);
                            }
                        },
                        NodeType::Solid => {
                            inner_no_slip(node, Q_EAST);
                        }
                    });
            }
            BoundaryFace::West => {
                self.nodes
                    .par_iter_mut()
                    .filter(|node| {
                        node.boundary_faces
                            .as_ref()
                            .unwrap()
                            .contains(&boundary_face)
                    })
                    .for_each(|node| match node.node_type {
                        NodeType::Fluid => match (density, velocity) {
                            (None, [Some(ux), Some(uy)]) => {
                                let f = &mut node.f;
                                let rho = 1.0 / (1.0 - ux)
                                    * (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7]));
                                f[1] = f[3] + (2.0 / 3.0) * rho * ux;
                                f[5] = f[7] - 0.5 * (f[2] - f[4])
                                    + 0.5 * rho * uy
                                    + (1.0 / 6.0) * rho * ux;
                                f[8] = f[6] + 0.5 * (f[2] - f[4]) - 0.5 * rho * uy
                                    + (1.0 / 6.0) * rho * ux;
                            }
                            (Some(rho), [None, Some(uy)]) => {
                                let f = &mut node.f;
                                let ux = 1.0
                                    - 1.0 / rho * (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7]));
                                f[1] = f[3] + (2.0 / 3.0) * rho * ux;
                                f[5] = f[7] - 0.5 * (f[2] - f[4])
                                    + 0.5 * rho * uy
                                    + (1.0 / 6.0) * rho * ux;
                                f[8] = f[6] + 0.5 * (f[2] - f[4]) - 0.5 * rho * uy
                                    + (1.0 / 6.0) * rho * ux;
                            }
                            _ => {
                                panic!("{}", ZOUHE_ERROR_MESSAGE);
                            }
                        },
                        NodeType::Solid => {
                            inner_no_slip(node, Q_WEST);
                        }
                    });
            }
            BoundaryFace::North => {
                self.nodes
                    .par_iter_mut()
                    .filter(|node| {
                        node.boundary_faces
                            .as_ref()
                            .unwrap()
                            .contains(&boundary_face)
                    })
                    .for_each(|node| match node.node_type {
                        NodeType::Fluid => match (density, velocity) {
                            (None, [Some(ux), Some(uy)]) => {
                                let f = &mut node.f;
                                let rho = 1.0 / (1.0 + uy)
                                    * (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6]));
                                f[4] = f[2] - (2.0 / 3.0) * rho * uy;
                                f[7] = f[5] + 0.5 * (f[1] - f[3])
                                    - 0.5 * rho * ux
                                    - (1.0 / 6.0) * rho * uy;
                                f[8] = f[6] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux
                                    - (1.0 / 6.0) * rho * uy;
                            }
                            (Some(rho), [Some(ux), None]) => {
                                let f = &mut node.f;
                                let uy = 1.0 / rho
                                    * (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6]))
                                    - 1.0;
                                f[4] = f[2] - (2.0 / 3.0) * rho * uy;
                                f[7] = f[5] + 0.5 * (f[1] - f[3])
                                    - 0.5 * rho * ux
                                    - (1.0 / 6.0) * rho * uy;
                                f[8] = f[6] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux
                                    - (1.0 / 6.0) * rho * uy;
                            }
                            _ => {
                                panic!("{}", ZOUHE_ERROR_MESSAGE);
                            }
                        },
                        NodeType::Solid => {
                            inner_no_slip(node, Q_NORTH);
                        }
                    });
            }
            BoundaryFace::South => {
                self.nodes
                    .par_iter_mut()
                    .filter(|node| {
                        node.boundary_faces
                            .as_ref()
                            .unwrap()
                            .contains(&boundary_face)
                    })
                    .for_each(|node| match node.node_type {
                        NodeType::Fluid => match (density, velocity) {
                            (None, [Some(ux), Some(uy)]) => {
                                let f = &mut node.f;
                                let rho = 1.0 / (1.0 - uy)
                                    * (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8]));
                                f[2] = f[4] + (2.0 / 3.0) * rho * uy;
                                f[5] = f[7] - 0.5 * (f[1] - f[3])
                                    + 0.5 * rho * ux
                                    + (1.0 / 6.0) * rho * uy;
                                f[6] = f[8] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux
                                    + (1.0 / 6.0) * rho * uy;
                            }
                            (Some(rho), [Some(ux), None]) => {
                                let f = &mut node.f;
                                let uy = 1.0
                                    - 1.0 / rho * (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8]));
                                f[2] = f[4] + (2.0 / 3.0) * rho * uy;
                                f[5] = f[7] - 0.5 * (f[1] - f[3])
                                    + 0.5 * rho * ux
                                    + (1.0 / 6.0) * rho * uy;
                                f[6] = f[8] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux
                                    + (1.0 / 6.0) * rho * uy;
                            }
                            _ => {
                                panic!("{}", ZOUHE_ERROR_MESSAGE);
                            }
                        },
                        NodeType::Solid => {
                            inner_no_slip(node, Q_SOUTH);
                        }
                    });
            }
        }
    }
}

fn inner_no_slip(node: &mut Node, q_face: [usize; 3]) {
    for q in q_face {
        node.f[Q_BAR[q]] = node.f_star[q];
    }
}
