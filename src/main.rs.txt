#![allow(unused_imports)]
use genetic_optimization::prelude::*;

fn main() {
    // polynomial a³ + b³ + c³ + j² + k² + l² + x + y + z
    let genome = Genome::new()
        .add_chromosome("cubes", Chromosome::new()
            .add_gene("a", Gene::new_with_range(1.0, -100.0, 100.0))
            .add_gene("b", Gene::new_with_range(1.0, -100.0, 100.0))
            .add_gene("c", Gene::new_with_range(1.0, -100.0, 100.0)))
        .add_chromosome("squares", Chromosome::new()
            .add_gene("j", Gene::new_with_range(1.0, -100.0, 100.0))
            .add_gene("k", Gene::new_with_range(1.0, -100.0, 100.0))
            .add_gene("l", Gene::new_with_range(1.0, -100.0, 100.0)))
        .add_chromosome("lines", Chromosome::new()
            .add_gene("x", Gene::new_with_range(1.0, -100.0, 100.0))
            .add_gene("y", Gene::new_with_range(1.0, -100.0, 100.0))
            .add_gene("z", Gene::new_with_range(1.0, -100.0, 100.0)))
        .build();

    Simulation::new()
        .genome(&genome)
        .eval(close_to_42)
        .print_settings(PrintSettings::PrintFull)
        .run(1000);
}

fn close_to_42(genome: &Genome) -> f64 {
    let a = genome.gene("cubes", "a").unwrap().value();
    let b = genome.gene("cubes", "b").unwrap().value();
    let c = genome.gene("cubes", "c").unwrap().value();
    let j = genome.gene("squares", "j").unwrap().value();
    let k = genome.gene("squares", "k").unwrap().value();
    let l = genome.gene("squares", "l").unwrap().value();
    let x = genome.gene("lines", "x").unwrap().value();
    let y = genome.gene("lines", "y").unwrap().value();
    let z = genome.gene("lines", "z").unwrap().value();

    let sum = a.powi(3) + b.powi(3) + c.powi(3) + j.powi(2) + k.powi(2) + l.powi(2) + x + y + z;

    // score is correlated to distance from 42
    // best possible score is 1.0
    -0.1 * (sum - 42.0).abs() + 1.0
}