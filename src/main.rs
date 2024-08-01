#![allow(unused_imports)]
use genetic_optimization::prelude::*;

macro_rules! quick_hello {
    ( $(h),* ) => {{
        let h = "hello";
        
        h
    }};
}

macro_rules! pos {
    ($($a:expr),*; $($b:expr),*) => {{
        let (mut sum1, mut sum2) = (0, 0);
        let h = "hello";
        $(
            sum1 += $a;
        )*
        $(
            sum2 += $b;
        )*
        h
    }};
}

fn main() {
    // std::env::set_var("RUST_BACKTRACE", "full");

    // polynomial a³ + b³ + c³ + j² + k² + l² + x + y + z
    // let genome = Genome::new(vec![
    //     ("cubes", "a", Gene::new_with_range(1.0, -100.0, 100.0)),
    //     ("cubes", "b", Gene::new_with_range(1.0, -100.0, 100.0)),
    //     ("cubes", "c", Gene::new_with_range(1.0, -100.0, 100.0)),
    //     ("squares", "j", Gene::new_with_range(1.0, -100.0, 100.0)),
    //     ("squares", "k", Gene::new_with_range(1.0, -100.0, 100.0)),
    //     ("squares", "l", Gene::new_with_range(1.0, -100.0, 100.0)),
    //     ("lines", "x", Gene::new_with_range(1.0, -100.0, 100.0)),
    //     ("lines", "y", Gene::new_with_range(1.0, -100.0, 100.0)),
    //     ("lines", "z", Gene::new_with_range(1.0, -100.0, 100.0))
    // ]);

    // let optimized = genome.simulate(100, 10000, close_to_42, SimParams::default(), true);

    // println!("Species: {}", optimized);

    let res = pos!(1,2,3; 4,5,6);

    let e = quick_hello!(h);
    println!("{h}");
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

