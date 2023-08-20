#![allow(unused_imports)]
use genetic_optimization::prelude::*;

fn main() {
    // polynomial a³ + b³ + c³ + j² + k² + l² + x + y + z
    let genome = Genome::new(vec![
        ("cubes", "a", Gene::new_with_range(1.0, -100.0, 100.0)),
        ("cubes", "b", Gene::new_with_range(1.0, -100.0, 100.0)),
        ("cubes", "c", Gene::new_with_range(1.0, -100.0, 100.0)),
        ("squares", "j", Gene::new_with_range(1.0, -100.0, 100.0)),
        ("squares", "k", Gene::new_with_range(1.0, -100.0, 100.0)),
        ("squares", "l", Gene::new_with_range(1.0, -100.0, 100.0)),
        ("lines", "x", Gene::new_with_range(1.0, -100.0, 100.0)),
        ("lines", "y", Gene::new_with_range(1.0, -100.0, 100.0)),
        ("lines", "z", Gene::new_with_range(1.0, -100.0, 100.0))
    ]);

    let time = std::time::Instant::now();
    genome.simulate(100000, 10000000, close_to_42, SimHyperParams::default(), true);
    println!("{}", time.elapsed().as_secs_f32());

    // let params = Genome::load_or_create("params.txt", 
    //     || genome.hyper_simulate(10, close_to_42, SimHyperParams::double_hyper(), true)).unwrap();
    
    // println!("{params}");

    // let mut score = 0;
    // for _ in 0..1000 {
    //     let optimized = genome.simulate(100, 10000, close_to_42, SimHyperParams::default(), false);
    //     let hyper_optimized = genome.simulate(100, 10000, close_to_42, SimHyperParams::from_genome(&params), false);

    //     if close_to_42(&hyper_optimized) > close_to_42(&optimized) {
    //         score += 1;
    //     }
    //     else {
    //         score -= 1;
    //     }
    // }

    // println!("{score}");

    // let optimized = genome.simulate(100, 10000, close_to_42, SimHyperParams::default(), true);
    // println!("Species: {}, Score: {}", optimized, close_to_42(&optimized));
    
    // let hyper_optimized: Genome = genome.simulate(100, 10000, close_to_42, SimHyperParams::from_genome(&params), true);
    // println!("Species: {}, Score: {}", hyper_optimized, close_to_42(&hyper_optimized));
}

fn close_to_42(genome: &Genome) -> f32 {
    if genome.gene("cubes", "a").is_none() {
        println!("{genome}");
    }

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