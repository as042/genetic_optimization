use rand::Rng;

use crate::prelude::*;

impl Genome {
    /// Simulates natural selection to optimize `self` for the given task.
    /// 
    /// The evaluator must be a function that takes a `Genome` with the same structure as `self` and returns a score (the greater, the better).
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// // polynomial a³ + b³ + c³ + j² + k² + l² + x + y + z
    /// let genome = Genome::new(vec![
    ///     ("cubes", "a", Gene::new_with_range(1.0, -100.0, 100.0)),
    ///     ("cubes", "b", Gene::new_with_range(1.0, -100.0, 100.0)),
    ///     ("cubes", "c", Gene::new_with_range(1.0, -100.0, 100.0)),
    ///     ("squares", "j", Gene::new_with_range(1.0, -100.0, 100.0)),
    ///     ("squares", "k", Gene::new_with_range(1.0, -100.0, 100.0)),
    ///     ("squares", "l", Gene::new_with_range(1.0, -100.0, 100.0)),
    ///     ("lines", "x", Gene::new_with_range(1.0, -100.0, 100.0)),
    ///     ("lines", "y", Gene::new_with_range(1.0, -100.0, 100.0)),
    ///     ("lines", "z", Gene::new_with_range(1.0, -100.0, 100.0)),
    /// ]);
    /// 
    /// // this will tell the simulator how to evaluate species
    /// fn close_to_42(genome: &Genome) -> f32 {
    ///     // variable bindings
    ///     let a = genome.gene("cubes", "a").unwrap().value();
    ///     let b = genome.gene("cubes", "b").unwrap().value();
    ///     let c = genome.gene("cubes", "c").unwrap().value();
    ///     let j = genome.gene("squares", "j").unwrap().value();
    ///     let k = genome.gene("squares", "k").unwrap().value();
    ///     let l = genome.gene("squares", "l").unwrap().value();
    ///     let x = genome.gene("lines", "x").unwrap().value();
    ///     let y = genome.gene("lines", "y").unwrap().value();
    ///     let z = genome.gene("lines", "z").unwrap().value();
    ///
    ///     // polynomial from earlier
    ///     let sum = a.powi(3) + b.powi(3) + c.powi(3) + j.powi(2) + k.powi(2) + l.powi(2) + x + y + z;
    ///
    ///     // score is correlated to distance from 42
    ///     // best possible score is 1.0
    ///     -0.1 * (sum - 42.0).abs() + 1.0
    /// }
    /// 
    /// // simulate 100 generations of optimization 
    /// let optimized = genome.simulate(100, close_to_42);
    /// assert_eq!(optimized.genes().len(), genome.genes().len());
    /// ```
    #[inline]
    pub fn simulate(&self, generations: usize, evaluator: fn(&Genome) -> f32) -> Genome {
        if generations == 0 { return self.clone(); }

        // Generation 1
        let mut species = gen_random_species(self, 15, true);
        species.push(self.clone());

        // Generations 2+
        for _ in 1..generations {
            let mut new_species = Vec::default();

            species.sort_unstable_by(|a, b| evaluator(b).partial_cmp(&evaluator(a)).unwrap());

            new_species.push(species[0].clone());
            for mutate in [false, true] {
                new_species.push(species[0].mate(&species[1], mutate));
                new_species.push(species[0].mate(&species[2], mutate));
                new_species.push(species[0].mate(&species[3], mutate));
                new_species.push(species[1].mate(&species[2], mutate));
                new_species.push(species[1].mate(&species[3], mutate));
                new_species.push(species[2].mate(&species[3], mutate));
            }

            new_species.extend(gen_random_species(self, 3, false));

            species = new_species;
        }

        species.sort_unstable_by(|a, b| evaluator(b).partial_cmp(&evaluator(a)).unwrap());
        
        species[0].clone()
    }

    #[inline]
    fn mate(&self, other: &Genome, mutate: bool) -> Genome {
        let mut rng = rand::thread_rng();

        let mut offspring = self.clone();
        for chromo in other.chromosomes() {
            if rng.gen_bool(0.5) {
                let _ = offspring.set_chromosome(chromo.0, chromo.1.clone());
            }
        }

        if mutate {
            for gene in offspring.clone().genes() {
                if rng.gen_bool(0.1) {
                    let _ = offspring.set_gene(gene.0, gene.1, gene.2.mutate());
                }
            }
        }

        offspring
    }
}

#[inline]
fn gen_random_species(template: &Genome, num: usize, mutate: bool) -> Vec<Genome> {
    let mut rng = rand::thread_rng();

    let mut species = Vec::default();
    for _ in 0..num {
        let mut spec = template.clone();
        for gene in template.genes() {
            if mutate {
                let _ = spec.set_gene(gene.0, gene.1, gene.2.mutate());
            }
            else {
                let _ = spec.set_gene(gene.0, gene.1, Gene::new_with_range(rng.gen_range(gene.2.min()..gene.2.max()), gene.2.min(), gene.2.max()));
            }            
        }

        species.push(spec)
    }

    species
}