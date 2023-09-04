use std::{collections::HashMap, fmt::Display, fs::{OpenOptions, self}, path::Path, io::Write};
use serde::{Serialize, Deserialize};

use crate::prelude::*;

/// A list of chromosomes that make up a specific specimen.
/// 
/// # Examples
/// 
/// ```
/// use genetic_optimization::prelude::*;
/// 
/// let animal = Genome::new()
///     .add_chromosome("eyes", Chromosome::new()
///         .add_gene("green", Gene::new(0.5))
///         .add_gene("vision_quality", Gene::new(0.9)))
///     .add_chromosome("behavior", Chromosome::new()
///         .add_gene("aggressiveness", Gene::new(0.2))
///         .add_gene("food_motivation", Gene::new(1.0)))
///     .build();
/// 
/// let alpha_specimen = Simulation::new()
///     .genome(&animal)
///     .eval(survivability)
///     .run(100);
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Genome {
    chromosomes: HashMap<String, Chromosome>
}

impl Genome {
    /// Creates a new `Genome` with the given genes.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let animal = Genome::new()
    ///     .add_chromosome("eyes", Chromosome::new()
    ///         .add_gene("green", Gene::new(0.5))
    ///         .add_gene("vision_quality", Gene::new(0.9)))
    ///     .add_chromosome("behavior", Chromosome::new()
    ///         .add_gene("aggressiveness", Gene::new(0.2))
    ///         .add_gene("food_motivation", Gene::new(1.0)))
    ///     .build();
    /// 
    /// assert!(animal.genes().len() > 0);
    /// ```
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a chromosome to `self` and returns `self`.
    #[inline]
    pub fn add_chromosome(&mut self, chromo_name: impl Into<String>, chromo: &Chromosome) -> &mut Self {
        self.insert_chromosome(chromo_name.into(), chromo.clone());
        self
    }

    /// Builds `self` and returns an owned value.
    #[inline]
    pub fn build(&mut self) -> Self {
        self.to_owned()
    }

    /// Returns the specified chromosome, if it exists.
    #[inline]
    pub fn chromosome(&self, chromo_name: impl Into<String>) -> Option<&Chromosome> {
        self.chromosomes.get(&chromo_name.into())
    }

    /// Returns the specified gene, if it exists.
    #[inline]
    pub fn gene(&self, chromo_name: impl Into<String>, gene_name: impl Into<String>) -> Option<&Gene> {
        self.chromosome(chromo_name)?.gene(gene_name)
    }

    /// Returns a reference to `self.chromosomes`.
    #[inline]
    pub fn chromosomes(&self) -> &HashMap<String, Chromosome> {
        &self.chromosomes
    }
    
    /// Returns a list of all genes, their labels, and their containers' labels.
    #[inline]
    pub fn genes(&self) -> Vec<(String, String, Gene)> {
        let mut genes = Vec::default();

        for chromo in self.chromosomes() {
            for gene in chromo.1.genes() {
                genes.push((chromo.0.to_owned(), gene.0.to_owned(), gene.1.clone()));
            }
            
        }        

        genes
    }

    /// Inserts a `Gene` into `self`.
    #[inline]
    pub fn insert_chromosome(&mut self, chromo_name: impl Into<String>, chromo: Chromosome) {
        self.chromosomes.insert(chromo_name.into(), chromo);
    }

    /// Inserts all given chromosomes into `self`.
    #[inline]
    pub fn insert_chromosomes(&mut self, chromos: Vec<(impl Into<String>, Chromosome)>) {
        for chromo in chromos {
            self.insert_chromosome(chromo.0, chromo.1);
        }
    }

    /// Inserts all given genes into `self`.
    #[inline]
    pub fn insert_genes(&mut self, genes: Vec<(impl Into<String> + Clone, impl Into<String>, Gene)>) {
        for gene in genes {
            if !self.chromosomes.contains_key(&gene.0.clone().into()) {
                self.insert_chromosome(&gene.0.clone().into(), Chromosome::new());
            }

            self.chromosomes.get_mut(&gene.0.into()).unwrap().insert_gene(&gene.1.into(), gene.2);
        }
    }

    /// Sets the value of the specified chromosome.
    #[inline]
    pub fn set_chromosome(&mut self, chromo_name: impl Into<String>, chromo: Chromosome) -> Result<(), String> {
        *self.chromosomes.get_mut(&chromo_name.into()).ok_or("No chromosome by that name exists")? = chromo;

        Ok(())
    }

    /// Sets the value of the specified gene.
    #[inline]
    pub fn set_gene(&mut self, chromo_name: impl Into<String>, gene_name: impl Into<String>, gene: Gene) -> Result<(), impl Into<String>> {
        self.chromosomes.get_mut(&chromo_name.into()).ok_or("No chromosome by that name exists")?.set_gene(gene_name, gene)
    }

    /// Saves `self` to file.
    #[inline]
    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), String> {
        let data = toml::to_string_pretty(self).ok().ok_or("Could not convert to toml")?;

        let mut file = OpenOptions::new()
            .create(true)
            .write(true)
            .open(path)
            .ok().ok_or("Could not open file")?;
        file.write(data.as_bytes()).ok().ok_or("Could not write to file")?;

        Ok(())
    }

    /// Loads file and returns its `Genome`.
    #[inline]
    pub fn load(path: impl AsRef<Path>) -> Result<Self, &'static str> {
        let string = fs::read_to_string(path).ok().ok_or("Could not open file")?;

        toml::from_str(&string).ok().ok_or("Invalid toml")
    }
    
    /// Attempts to load `Genome`, if fails, uses existing.
    #[inline]
    pub fn load_or(path: impl AsRef<Path>, genome: Self) -> Self {
        if let Ok(res) = Genome::load(path) {
            return res;
        }
        
        genome
    }

    /// Attempts to load `Genome`, if fails, creates a new `Genome`.
    #[inline]
    pub fn load_or_else(path: impl AsRef<Path>, genome: impl FnOnce() -> Self) -> Self {
        if let Ok(res) = Genome::load(path) {
            return res;
        }
        
        genome()
    }

    /// Attempts to load `Genome`, if fails, creates a new `Genome` and saves it to the file.
    #[inline]
    pub fn load_or_create(path: impl AsRef<Path> + Clone, genome: impl FnOnce() -> Self) -> Result<Self, &'static str> {
        if let Ok(res) = Genome::load(path.clone()) {
            return Ok(res);
        }
        
        let genome = genome();
        genome.save(path).ok().ok_or("Failed to save")?;

        Ok(genome)
    }
}

impl Display for Genome {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut display_obj: HashMap<String, HashMap<String, f64>> = HashMap::new();
        for chromo in &self.chromosomes {
            display_obj.insert(chromo.0.to_string(), HashMap::new());
            for gene in chromo.1.genes() {
                display_obj.get_mut(chromo.0).unwrap().insert(gene.0.to_string(), gene.1.value());
            }
        }
        write!(f, "{:#?}", display_obj)
    }
}