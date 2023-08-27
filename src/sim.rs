use std::{sync::mpsc::{Sender, Receiver, self}, thread};

use rand::Rng;
use threadpool::ThreadPool;

use crate::prelude::*;

/// Hyper parameters for `simulate()`.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct SimHyperParams {
    pub elitism_survivors: usize,
    pub elitism_reproducers: usize,
    pub random_reproducers: usize,
    pub num_random_species: usize,

    pub crossover_chance_per_gene: f32,
    pub offspring_mutation_chance: f32,
    pub offspring_mutation_randomness_weight: f32,
    pub random_species_randomness_weight: f32,
}

impl SimHyperParams {
    #[allow(dead_code)]
    #[deprecated]
    #[inline]
    fn debug() -> Self {
        SimHyperParams {
            elitism_survivors: 1,
            elitism_reproducers: 0,
            random_reproducers: 0,
            num_random_species: 1,
            crossover_chance_per_gene: 0.5,
            offspring_mutation_chance: 0.5,
            offspring_mutation_randomness_weight: 0.5,
            random_species_randomness_weight: 0.5,
        }
    }

    /// Returns the default `SimHyperParams` for use in `hyper_simulate()`.
    #[inline]
    pub fn double_hyper() -> Self {
        SimHyperParams {
            elitism_survivors: 3,
            elitism_reproducers: 3,
            random_reproducers: 2,
            num_random_species: 2,
            crossover_chance_per_gene: 0.1,
            offspring_mutation_chance: 0.1,
            offspring_mutation_randomness_weight: 0.25,
            random_species_randomness_weight: 0.5,
        }
    }

    /// Creates a `SimHyperParams` from `Genome`.
    #[inline]
    pub fn from_genome(genome: &Genome) -> Self {
        let elitism_survivors = genome.gene("species", "elitism_survivors").unwrap().value();
        let elitism_reproducers = genome.gene("species", "elitism_reproducers").unwrap().value();
        let random_reproducers = genome.gene("species", "random_reproducers").unwrap().value();
        let num_random_species = genome.gene("species", "num_random_species").unwrap().value();
        let crossover_chance_per_gene = genome.gene("recombination", "crossover_chance_per_gene").unwrap().value();
        let offspring_mutation_chance = genome.gene("recombination", "offspring_mutation_chance").unwrap().value();
        let offspring_mutation_randomness_weight = genome.gene("recombination", "offspring_mutation_randomness_weight").unwrap().value();
        let random_species_randomness_weight = genome.gene("random_species", "random_species_randomness_weight").unwrap().value();
    
        SimHyperParams {
            elitism_survivors: elitism_survivors.round() as usize,
            elitism_reproducers: elitism_reproducers.round() as usize,
            random_reproducers: random_reproducers.round() as usize,
            num_random_species: num_random_species.round() as usize,
            crossover_chance_per_gene,
            offspring_mutation_chance,
            offspring_mutation_randomness_weight,
            random_species_randomness_weight,
        }
    }

    /// Returns the number of species a generation will have if a simulation is run using `self`.
    #[inline]
    pub fn species_per_generation(&self) -> usize {
        self.elitism_survivors + self.num_random_species + self.random_reproducers + (self.elitism_reproducers as i32 * (self.elitism_reproducers as i32 - 1) / 2) as usize
    }
}

impl Default for SimHyperParams {
    #[inline]
    fn default() -> Self {
        Self { 
            elitism_survivors: 1, 
            elitism_reproducers: 30, 
            random_reproducers: 11,
            num_random_species: 10, 
            crossover_chance_per_gene: 0.1, 
            offspring_mutation_chance: 0.1, 
            offspring_mutation_randomness_weight: 0.25, 
            random_species_randomness_weight: 1.0
        }
    }
}

impl Genome {
    /// Trains the GA itself to best be able to tackle the given problem.
    #[inline]
    pub fn hyper_simulate(&self, generations: usize, evaluator: fn(&Genome) -> f32, hyper_params: SimHyperParams, print: bool) -> Genome {
        let hyper_genome = Genome::new(vec![
            ("species", "elitism_survivors", Gene::new_with_range(1.0, 0.0, 100.0)),
            ("species", "elitism_reproducers", Gene::new_with_range(13.0, 0.0, 20.0)),
            ("species", "random_reproducers", Gene::new_with_range(11.0, 0.0, 100.0)),
            ("species", "num_random_species", Gene::new_with_range(10.0, 0.0, 100.0)),
            ("recombination", "crossover_chance_per_gene", Gene::new_with_range(0.1, 0.0, 1.0)),
            ("recombination", "offspring_mutation_chance", Gene::new_with_range(0.1, 0.0, 1.0)),
            ("recombination", "offspring_mutation_randomness_weight", Gene::new_with_range(0.25, 0.0, 1.0)),
            ("random_species", "random_species_randomness_weight", Gene::new_with_range(1.0, 0.0, 1.0))
        ]);

        hyper_genome.hyper_sim(self, generations, 1000000, evaluator, hyper_params, print)
    }

    #[inline]
    fn hyper_sim(&self, template: &Genome, mut generations: usize, species_limit: usize, evaluator: fn(&Genome) -> f32, hyper_params: SimHyperParams, print: bool) -> Genome {
        if generations == 0 { return self.clone(); }

        let mut rng = rand::thread_rng();

        // Generation 1
        let mut species = gen_random_species(self, hyper_params.species_per_generation() - 1, hyper_params.random_species_randomness_weight);
        species.push(self.clone());

        // Generations 2+
        for i in 1..generations {
            let mut new_species = Vec::default();

            species.sort_unstable_by(|a, b| hyper_eval(b, template, evaluator).partial_cmp(&hyper_eval(a, template, evaluator)).unwrap());

            if print { println!("Generation: {i}, Score: {}", hyper_eval(&species[0], template, evaluator)); }

            // Elitism
            for j in 0..hyper_params.species_per_generation() {
                // Elitism survivors
                if j < hyper_params.elitism_survivors {
                    new_species.push(species[j].clone());
                }
                // Elitism reproducers
                if j + 1 < hyper_params.elitism_reproducers {
                    for k in j + 1..hyper_params.elitism_reproducers {
                        new_species.push(species[j].mate(&species[k], hyper_params.crossover_chance_per_gene, hyper_params.offspring_mutation_chance, hyper_params.offspring_mutation_randomness_weight));
                    }
                }
                if j >= hyper_params.elitism_survivors && j + 1 >= hyper_params.elitism_reproducers {
                    break;
                }
            }

            // Random reproducers
            for _ in 0..hyper_params.random_reproducers {
                let idx = rng.gen_range(1..hyper_params.species_per_generation());
                new_species.push(species[0].mate(&species[idx], hyper_params.crossover_chance_per_gene, hyper_params.offspring_mutation_chance, hyper_params.offspring_mutation_randomness_weight));
            }

            new_species.extend(gen_random_species(self, hyper_params.num_random_species, hyper_params.random_species_randomness_weight));

            species = new_species;

            generations = i;

            // Species limit
            if (i + 1) * hyper_params.species_per_generation() >= species_limit {
                break;
            }
        }

        species.sort_unstable_by(|a, b| hyper_eval(b, template, evaluator).partial_cmp(&hyper_eval(a, template, evaluator)).unwrap());

        if print { println!("Generation: {}, Score: {}", generations + 1, hyper_eval(&species[0], template, evaluator)); }
        
        species[0].clone()
    }

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
    pub fn simulate(&self, mut generations: usize, species_limit: usize, evaluator: fn(&Genome) -> f32, hyper_params: SimHyperParams, print: bool) -> Genome {
        if generations == 0 { return self.clone(); }

        let mut rng = rand::thread_rng();

        // Generation 1
        let mut species = gen_random_species(self, hyper_params.species_per_generation() - 1, hyper_params.random_species_randomness_weight);
        species.push(self.clone());

        // Generations 2+
        for i in 1..generations {
            let mut new_species = Vec::default();

            // sort_species(&mut species, evaluator);
            species.sort_by_cached_key(|x| (evaluator(x) * -1_000_000.0) as i64);

            if print { println!("Generation: {i}, Score: {}", evaluator(&species[0])); }

            // Elitism
            for j in 0..hyper_params.species_per_generation() {
                // Elitism survivors
                if j < hyper_params.elitism_survivors {
                    new_species.push(species[j].clone());
                }
                // Elitism reproducers
                if j + 1 < hyper_params.elitism_reproducers {
                    for k in j + 1..hyper_params.elitism_reproducers {
                        new_species.push(species[j].mate(&species[k], hyper_params.crossover_chance_per_gene, hyper_params.offspring_mutation_chance, hyper_params.offspring_mutation_randomness_weight));
                    }
                }
                if j >= hyper_params.elitism_survivors && j + 1 >= hyper_params.elitism_reproducers {
                    break;
                }
            }

            // Random reproducers
            for _ in 0..hyper_params.random_reproducers {
                let idx = rng.gen_range(1..hyper_params.species_per_generation());
                new_species.push(species[0].mate(&species[idx], hyper_params.crossover_chance_per_gene, hyper_params.offspring_mutation_chance, hyper_params.offspring_mutation_randomness_weight));
            }

            new_species.extend(gen_random_species(self, hyper_params.num_random_species, hyper_params.random_species_randomness_weight));

            species = new_species;

            generations = i;

            // Species limit
            if (i + 1) * hyper_params.species_per_generation() >= species_limit {
                break;
            }
        }

        // species.sort_by_cached_key(|x| (evaluator(x) * -1_000_000.0) as i64);
        sort_species(&mut species, evaluator);

        if print { println!("Generation: {}, Score: {}", generations + 1, evaluator(&species[0])); }
        
        species[0].clone()
    }

    #[inline]
    fn mate(&self, other: &Genome, crossover_chance: f32, mutation_chance: f32, randomness_weight: f32) -> Genome {
        let mut rng = rand::thread_rng();

        let mut offspring = self.clone();

        // Recombination
        for chromo in offspring.clone().chromosomes() {
            let mut use_other = rng.gen_bool(0.5);
            for g in chromo.1.genes() {
                // Pick which chromosome to take gene from
                let mut gene = g.1;
                if use_other {
                    gene = other.gene(chromo.0, g.0).unwrap();
                }
    
                // Mutation
                if rng.gen_bool(mutation_chance as f64) {
                    let _ = offspring.set_gene(chromo.0, g.0, gene.mutate(randomness_weight));
                }
    
                // Crossover
                if rng.gen_bool(crossover_chance as f64) {
                    use_other = !use_other;
                }
            }
        }


        offspring
    }
}

#[inline]
fn sort_species(species: &mut Vec<Genome>, eval: fn(&Genome) -> f32) {
    let pool = ThreadPool::new(16);

    let (tx, rx)= mpsc::channel();

    for k in 0..species.len() {
        let thread_tx = tx.clone();
        let spec = species[k].clone();

        pool.execute(move || {
            thread_tx.send((k, eval(&spec))).unwrap();
        });
    }

    let mut new = Vec::with_capacity(species.len());
    for _ in 0..species.len() {
        let recv = rx.recv().unwrap();
        new.push((species[recv.0].clone(), recv.1));
    }
    
    // pool.join();

    new.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    *species = new.iter().map(|x| x.0.clone()).collect();
}

#[inline]
fn gen_random_species(template: &Genome, num: usize, randomness_weight: f32) -> Vec<Genome> {
    let mut species = Vec::default();
    for _ in 0..num {
        let mut spec = template.clone();
        for gene in template.genes() {
            let _ = spec.set_gene(gene.0, gene.1, gene.2.mutate(randomness_weight));         
        }

        species.push(spec)
    }

    species
}

#[inline]
fn hyper_eval(genome: &Genome, template: &Genome, eval: fn(&Genome) -> f32) -> f32 {
    let hyper_params = SimHyperParams::from_genome(genome);

    let mut score = 0.0;
    for _ in 0..32 {
        score += eval(&template.simulate(50, 10000, eval, hyper_params, false));
    }

    score / 32.0
}