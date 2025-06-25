use clap::{arg, command, value_parser, Command};
use lbm_lesc as lbm;
use rayon::ThreadPoolBuilder;
use std::env;

fn main() {
    let matches = command!()
        .arg(
            arg!(
                -d --dimensions <DIMENSIONS> "Sets the dimension of the simulation: 2 or 3"
            )
            .required(true)
            .value_parser(value_parser!(usize)),
        )
        .arg(
            arg!(
                -n --number_of_threads <NUMBER_OF_THREADS> "Sets the number of threads: 1, 2, 4, 8, 16 or 32"
            )
            .required(true)
            .value_parser(value_parser!(usize)),
        )
        .subcommand(
            Command::new("run")
                .about("Runs the simulation")
                .arg(
                    arg!(
                        -b --benchmark "Runs the benchmark"
                    )
                    .required(false),
                ),
        )
        .subcommand(
            Command::new("post")
                .about("Runs the post-processing: writes the vtk files")
        )
        .get_matches();

    if let Some(&num_threads) = matches.get_one::<usize>("number_of_threads") {
        ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .unwrap();
    }

    match matches.subcommand() {
        Some(("run", sub_matches)) => {
            if let Some(dimensions) = matches.get_one::<usize>("dimensions") {
                match dimensions {
                    2 => match sub_matches.get_one::<bool>("benchmark").unwrap() {
                        false => lbm::d2q9::run(),
                        true => lbm::d2q9::run_benchmark(),
                    },
                    3 => match sub_matches.get_one::<bool>("benchmark").unwrap() {
                        false => lbm::d3q27::run(),
                        true => lbm::d3q27::run_benchmark(),
                    },
                    _ => {
                        eprintln!(
                            "Error: the number of dimensions {dimensions} is not valid. Please, use 2 or 3."
                        );
                        std::process::exit(1);
                    }
                }
            }
        }
        Some(("post", _)) => {
            if let Some(dimensions) = matches.get_one::<usize>("dimensions") {
                match dimensions {
                    2 => lbm::d2q9::post::vtk::run_vtk_post_processing().unwrap(),
                    3 => lbm::d3q27::post::vtk::run_vtk_post_processing().unwrap(),
                    _ => {
                        eprintln!(
                            "Error: the number of dimensions {dimensions} is not valid. Please, use 2 or 3."
                        );
                        std::process::exit(1);
                    }
                }
            }
        }
        _ => {}
    }
}
