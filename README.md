# genetic_optimization

Genetic algorithm framework for solving optimization problems. Parameters are described as a
`Genome` built from named `Chromosome`s of named `Gene`s, where each gene is a single `f64` with an
optional min and max. You supply an evaluation function that scores a genome (higher is better), and
the simulation evolves a population toward it: each generation sorts the species by score, keeps the
elites, mates them with each other and with random partners using per-gene crossover and mutation,
and refills the rest of the population with fresh random species. Mutation rates, elitism counts,
and population size are all exposed through `SimHyperParams`. Evaluation can run single-threaded or
across threads, and the `Auto` setting picks between them by timing both on a sample generation.
Genomes serialize to TOML, so an optimized result can be saved and reloaded with `save`, `load`, and
`load_or_create`.

## Example

```rust
use genetic_optimization::prelude::*;

fn main() {
    let genome = Genome::new()
        .add_chromosome("cubes", Chromosome::new()
            .add_gene("a", Gene::new_with_range(1.0, -100.0, 100.0))
            .add_gene("b", Gene::new_with_range(1.0, -100.0, 100.0)))
        .add_chromosome("lines", Chromosome::new()
            .add_gene("x", Gene::new_with_range(1.0, -100.0, 100.0)))
        .build();

    let optimized = Simulation::new()
        .genome(&genome)
        .eval(close_to_42)
        .print_settings(PrintSettings::PrintScores)
        .run(100);

    println!("{}", optimized);
    optimized.save("best.toml").unwrap();
}

// a³ + b³ + x, scored by how close the result is to 42. Best possible score is 1.0.
fn close_to_42(genome: &Genome) -> f64 {
    let a = genome.gene("cubes", "a").unwrap().value();
    let b = genome.gene("cubes", "b").unwrap().value();
    let x = genome.gene("lines", "x").unwrap().value();

    let sum = a.powi(3) + b.powi(3) + x;

    -0.1 * (sum - 42.0).abs() + 1.0
}
```