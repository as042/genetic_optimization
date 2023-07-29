#![allow(unused_imports)]
use genetic_optimization::prelude::*;

fn main() {
    let genome = Genome::new(vec![
        ("chromo1", "gene1", Gene::new(1.0))
    ]);

    println!("{:#?}", genome);
}