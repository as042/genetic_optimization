use std::{sync::mpsc, thread, time::Instant};
use rand::Rng;

use crate::prelude::*;

#[derive(Clone, Default)]
pub struct Simulation {
    genome: Option<Genome>,
    species_limit: usize,
    eval: Option<fn(&Genome) -> f64>,
    eval_with_util: Option<fn(&Genome, fn(Vec<f64>) -> f64) -> f64>,
    util: Option<fn(Vec<f64>) -> f64>,
    hyper_params: SimHyperParams,
    parallelism: Parallelism,
    print_settings: PrintSettings,
}

impl Simulation {
    /// Creates a new `Simulation` that can be configured.
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a genome to `self`.
    #[inline]
    pub fn genome(&mut self, genome: &Genome) -> &mut Self {
        self.genome = Some(genome.clone());
        self
    }

    /// Adds a species limit to `self`.
    #[inline]
    pub fn species_limit(&mut self, species_limit: usize) -> &mut Self {
        self.species_limit = species_limit;
        self
    }

    /// Adds an evaluation function to `self`. An evaluation function is required.
    /// 
    /// # Panics
    /// Will panic if `self.eval_with_util()` has already been run.
    #[inline]
    pub fn eval(&mut self, eval: fn(&Genome) -> f64) -> &mut Self {
        assert!(self.eval_with_util.is_none());

        self.eval = Some(eval);
        self
    }

    /// Adds an evaluation function with an additional util function to `self`. An evaluation function is required.
    /// 
    /// # Panics
    /// Will panic if `self.eval()` has already been run.
    #[inline]
    pub fn eval_with_util(&mut self, eval: fn(&Genome, fn(Vec<f64>) -> f64) -> f64, util: fn(Vec<f64>) -> f64) -> &mut Self {
        assert!(self.eval.is_none());

        self.eval_with_util = Some(eval);
        self.util = Some(util);
        self
    }

    /// Sets the hyper parameters of `self`.
    #[inline]
    pub fn hyper_params(&mut self, hyper_params: SimHyperParams) -> &mut Self {
        self.hyper_params = hyper_params;
        self
    }

    /// Sets the parallelism option of `self`.
    #[inline]
    pub fn parallelism(&mut self, parallelism: Parallelism) -> &mut Self {
        self.parallelism = parallelism;
        self
    }

    /// Sets the print settings of `self`.
    #[inline]
    pub fn print_settings(&mut self, print_settings: PrintSettings) -> &mut Self {
        self.print_settings = print_settings;
        self
    }

    /// Runs the simulation using the configuration in `self`, optimizing the parameter for the specified task.
    /// The evaluator must be a function that takes a `Genome` with the same structure as `self` and returns a score (the greater, the better).
    /// 
    /// # Panics
    /// 
    /// This function will panic if no genome or evaluation function is given.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// // polynomial a³ + b³ + c³ + j² + k² + l² + x + y + z
    /// let genome = Genome::new()
    ///     .add_chromosome("cubes", Chromosome::new()
    ///         .add_gene("a", Gene::new_with_range(1.0, -100.0, 100.0))
    ///         .add_gene("b", Gene::new_with_range(1.0, -100.0, 100.0))
    ///         .add_gene("c", Gene::new_with_range(1.0, -100.0, 100.0)))
    ///     .add_chromosome("squares", Chromosome::new()
    ///         .add_gene("j", Gene::new_with_range(1.0, -100.0, 100.0))
    ///         .add_gene("k", Gene::new_with_range(1.0, -100.0, 100.0))
    ///         .add_gene("l", Gene::new_with_range(1.0, -100.0, 100.0)))
    ///     .add_chromosome("lines", Chromosome::new()
    ///         .add_gene("x", Gene::new_with_range(1.0, -100.0, 100.0))
    ///         .add_gene("y", Gene::new_with_range(1.0, -100.0, 100.0))
    ///         .add_gene("z", Gene::new_with_range(1.0, -100.0, 100.0)))
    ///     .build();
    /// 
    /// // this will tell the simulator how to evaluate species
    /// fn close_to_42(genome: &Genome) -> f64 {
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
    /// let optimized = Simulation::new()
    ///     .genome(&genome)
    ///     .eval(close_to_42)
    ///     .print_settings(PrintSettings::PrintFull)
    ///     .run(100);
    /// 
    /// assert_eq!(optimized.genes().len(), genome.genes().len());
    /// ```
    #[inline]
    pub fn run(&self, generations: usize) -> Genome {
        assert!(self.genome.is_some());
        assert_ne!(self.eval.is_some(), self.eval_with_util.is_some());

        self.genome.clone().unwrap().sim(generations, self.species_limit, self.eval, self.eval_with_util, self.util, self.hyper_params, self.parallelism, self.print_settings)
    }
}

/// Hyper parameters for `simulate()`.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct SimHyperParams {
    pub elitism_survivors: usize,
    pub elitism_reproducers: usize,
    pub random_reproducers: usize,
    pub num_random_species: usize,

    pub crossover_chance_per_gene: f64,
    pub offspring_mutation_chance: f64,
    pub offspring_mutation_randomness_weight: f64,
    pub random_species_randomness_weight: f64,
}

impl SimHyperParams {
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
            elitism_reproducers: 13, 
            random_reproducers: 11,
            num_random_species: 10, 
            crossover_chance_per_gene: 0.1, 
            offspring_mutation_chance: 0.1, 
            offspring_mutation_randomness_weight: 0.25, 
            random_species_randomness_weight: 1.0
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub enum Parallelism {
    Single,
    Multi,
    #[default]
    Auto,
}

impl Parallelism {
    #[inline]
    fn is_multi(
        &self, 
        eval: Option<fn(&Genome) -> f64>, 
        eval_with_util: Option<fn(&Genome, fn(Vec<f64>) -> f64) -> f64>, 
        util: Option<fn(Vec<f64>) -> f64>, 
        template: &Genome, 
        species_per_generation: usize
    ) -> bool 
    {
        if self == &Self::Single { return false; }
        else if self == &Self::Multi { return true; }
        else {
            let mut species = gen_random_species(template, species_per_generation, 0.5);
            let mut species2 = species.clone();

            let inst = Instant::now();

            if let Some(eval) = eval {
                species.sort_by_cached_key(|x| (eval(x) * -1_000_000.0) as i64);
            }
            else {
                species.sort_by_cached_key(|x| (eval_with_util.unwrap()(x, util.unwrap()) * -1_000_000.0) as i64);
            }

            let time = inst.elapsed().as_nanos();
            let inst = Instant::now();

            if let Some(eval) = eval {
                sort_species(&mut species2, eval);
            }
            else {
                sort_species_util(&mut species2, eval_with_util.unwrap(), util.unwrap());
            }

            let time2 = inst.elapsed().as_nanos();

            if time > time2 {
                return true;
            }
            else {
                return false;
            }
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub enum PrintSettings {
    #[default]
    NoPrint,
    PrintScores,
    PrintFull,
}

impl PrintSettings {
    #[inline]
    fn print_scores(&self) -> bool {
        if self == &Self::PrintScores || self == &Self::PrintFull {
            return true;
        }

        false
    }
}

impl Genome {
    /// Simulates natural selection to optimize `self` for the given task.
    /// 
    /// The evaluator must be a function that takes a `Genome` with the same structure as `self` and returns a score (the greater, the better).
    #[deprecated = "Use `Simulation::new()` instead"]
    #[inline]
    pub fn simulate(
        &self, 
        generations: usize, 
        species_limit: usize, 
        eval: Option<fn(&Genome) -> f64>, 
        eval_with_util: Option<fn(&Genome, fn(Vec<f64>) -> f64) -> f64>,
        util: Option<fn(Vec<f64>) -> f64>,
        hyper_params: SimHyperParams, 
        parallellism: Parallelism, 
        print: PrintSettings
    ) -> Genome 
    {
        self.sim(generations, species_limit, eval, eval_with_util, util, hyper_params, parallellism, print)
    }

    #[inline]
    pub fn sim(
        &self, 
        mut generations: usize, 
        mut species_limit: usize, 
        eval: Option<fn(&Genome) -> f64>, 
        eval_with_util: Option<fn(&Genome, fn(Vec<f64>) -> f64) -> f64>,
        util: Option<fn(Vec<f64>) -> f64>,
        hyper_params: SimHyperParams, 
        parallellism: Parallelism, 
        print: PrintSettings
    ) -> Genome 
    {
        assert_ne!(eval.is_some(), eval_with_util.is_some());
        assert_eq!(eval_with_util.is_some(), util.is_some());

        if generations == 0 { return self.clone(); }
        if species_limit == 0 { species_limit = usize::MAX }

        let mut rng = rand::thread_rng();

        // Parallelism
        let inst = Instant::now();
        if print == PrintSettings::PrintFull { println!("Determining parallelism option..."); }
        let multi = parallellism.is_multi(eval, eval_with_util, util, self, hyper_params.species_per_generation());
        if print == PrintSettings::PrintFull { 
            if multi {
                println!("Option `Multi` chosen after {}s", inst.elapsed().as_secs_f32());
            }
            else {
                println!("Option `Single` chosen after {}s", inst.elapsed().as_secs_f32());
            }
        }

        // Generation 1
        let mut species = gen_random_species(self, hyper_params.species_per_generation() - 1, hyper_params.random_species_randomness_weight);
        species.push(self.clone());

        // Generations 2+
        for i in 1..generations {
            let mut new_species = Vec::default();

            // Sorting
            if let Some(eval) = eval {
                if multi { sort_species(&mut species, eval); }
                else { species.sort_by_cached_key(|x| (eval(x) * -1_000_000.0) as i64); }
            }
            else {
                if multi { sort_species_util(&mut species, eval_with_util.unwrap(), util.unwrap()); }
                else { species.sort_by_cached_key(|x| (eval_with_util.unwrap()(x, util.unwrap()) * -1_000_000.0) as i64); }
            }

            // Printing
            if let Some(eval) = eval { if print.print_scores() { println!("Generation: {i}, Score: {}", eval(&species[0])); }}
            else { if print.print_scores() { println!("Generation: {i}, Score: {}", eval_with_util.unwrap()(&species[0], util.unwrap())); }}

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

        // Sorting
        if let Some(eval) = eval {
            if multi { sort_species(&mut species, eval); }
            else { species.sort_by_cached_key(|x| (eval(x) * -1_000_000.0) as i64); }
        }
        else {
            if multi { sort_species_util(&mut species, eval_with_util.unwrap(), util.unwrap()); }
            else { species.sort_by_cached_key(|x| (eval_with_util.unwrap()(x, util.unwrap()) * -1_000_000.0) as i64); }
        }

        // Printing
        if let Some(eval) = eval { if print.print_scores() { println!("Generation: {}, Score: {}", generations + 1, eval(&species[0])); }}
        else { if print.print_scores() { println!("Generation: {}, Score: {}", generations + 1, eval_with_util.unwrap()(&species[0], util.unwrap())); }}
        
        species[0].clone()
    }

    #[inline]
    fn mate(&self, other: &Genome, crossover_chance: f64, mutation_chance: f64, randomness_weight: f64) -> Genome {
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
fn sort_species(species: &mut Vec<Genome>, eval: fn(&Genome) -> f64) {
    let (tx, rx)= mpsc::channel();

    for k in 0..species.len() {
        let thread_tx = tx.clone();
        let spec = species[k].clone();

        thread::spawn(move || {
            thread_tx.send((k, eval(&spec))).unwrap();
        });
    }

    let mut new = Vec::with_capacity(species.len());
    for _ in 0..species.len() {
        let recv = rx.recv().unwrap();
        new.push((species[recv.0].clone(), recv.1));
    }

    new.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    *species = new.iter().map(|x| x.0.clone()).collect();
}

#[inline]
fn sort_species_util(species: &mut Vec<Genome>, eval: fn(&Genome, fn(Vec<f64>) -> f64) -> f64, util: fn(Vec<f64>) -> f64) {
    let (tx, rx)= mpsc::channel();

    for k in 0..species.len() {
        let thread_tx = tx.clone();
        let spec = species[k].clone();

        thread::spawn(move || {
            thread_tx.send((k, eval(&spec, util))).unwrap();
        });
    }

    let mut new = Vec::with_capacity(species.len());
    for _ in 0..species.len() {
        let recv = rx.recv().unwrap();
        new.push((species[recv.0].clone(), recv.1));
    }

    new.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    *species = new.iter().map(|x| x.0.clone()).collect();
}

#[inline]
fn gen_random_species(template: &Genome, num: usize, randomness_weight: f64) -> Vec<Genome> {
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