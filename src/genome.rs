use std::{collections::HashMap, fmt::Display, fs::{OpenOptions, self}, path::Path, io::Write};
use serde::{Serialize, Deserialize};

use crate::prelude::*;

/// A list of chromosomes that make up a specific specimen.
/// 
/// # Examples
/// 
/// ```
/// let animal = Genome::new(vec![
///     ("eyes", "green", Gene::new(0.5)),
///     ("eyes", "vision_quality", Gene::new(0.9)),
///     ("behavior", "aggressiveness", Gene::new(0.2)),
///     ("behavior", "food_motivation", Gene::new(1.0))
/// ]);
/// 
/// let alpha_specimen = animal.simulate(100, survivability_evaluator);
/// ```
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Genome {
    chromosomes: HashMap<String, Chromosome>,
}

impl Genome {
    /// Creates a new `Genome` with the given genes.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let genome = Genome::new(vec![("chromosome_label", "gene_label", Gene::default())]);
    /// ```
    #[inline]
    pub fn new(genes: Vec<(&str, &str, Gene)>) -> Self {
        let mut genome = Self { 
            ..Default::default()
        };

        genome.insert_genes(genes);

        genome
    }

    /// Returns the specified chromosome, if it exists.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let genome = Genome::new(vec![("chromosome_label", "gene_label", Gene::default())]);
    /// let chromo = genome.chromosome("chromosome_label");
    /// assert!(chromo.is_some());
    /// ```
    #[inline]
    pub fn chromosome(&self, chromo_name: &str) -> Option<&Chromosome> {
        self.chromosomes.get(chromo_name)
    }

    /// Returns the specified gene, if it exists.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let genome = Genome::new(vec![("chromosome_label", "gene_label", Gene::default())]);
    /// let gene = genome.gene("chromosome_label", "gene_label");
    /// assert!(gene.is_some());
    /// ```
    #[inline]
    pub fn gene(&self, chromo_name: &str, gene_name: &str) -> Option<&Gene> {
        self.chromosome(chromo_name)?.gene(gene_name)
    }

    /// Returns a reference to `self.chromosomes`.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let genome = Genome::new(vec![("chromosome_label", "gene_label", Gene::default())]);
    /// let chromosomes = genome.chromosomes();
    /// assert!(chromosomes.contains_key("chromosome_label"));
    /// ```
    #[inline]
    pub fn chromosomes(&self) -> &HashMap<String, Chromosome> {
        &self.chromosomes
    }
    
    /// Returns a list of all genes, their labels, and their containers' labels.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let genome = Genome::new(vec![("chromosome_label", "gene_label", Gene::default())]);
    /// let genes = genome.genes();
    /// assert_eq!(genes[0].2.value(), 0.0);
    /// ```
    #[inline]
    pub fn genes(&self) -> Vec<(&str, &str, Gene)> {
        let mut genes = Vec::default();

        for chromo in self.chromosomes() {
            for gene in chromo.1.genes() {
                genes.push((chromo.0.as_str(), gene.0.as_str(), gene.1.clone()));
            }
            
        }        

        genes
    }

    /// Inserts a `Gene` into `self`.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let mut genome = Genome::new(vec![("examples", "example1", Gene::default())]);
    /// genome.insert_chromosome("examples2", Chromosome::default());
    /// ```
    #[inline]
    pub fn insert_chromosome(&mut self, chromo_name: &str, chromo: Chromosome) {
        self.chromosomes.insert(chromo_name.to_string(), chromo);
    }

    /// Inserts all given chromosomes into `self`.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let mut genome = Genome::new(vec![("examples", "example1", Gene::default())]);
    /// genome.insert_chromosomes(vec![("examples2", Chromosome::default()), ("examples3", Chromosome::default())]);
    /// ```
    #[inline]
    pub fn insert_chromosomes(&mut self, chromos: Vec<(&str, Chromosome)>) {
        for chromo in chromos {
            self.insert_chromosome(chromo.0, chromo.1);
        }
    }

    /// Inserts all given genes into `self`.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let mut genome = Genome::default()
    /// genome.insert_genes(vec!["chromosome_label", "gene_label", Gene::default()]);
    /// ```
    #[inline]
    pub fn insert_genes(&mut self, genes: Vec<(&str, &str, Gene)>) {
        for gene in genes {
            if !self.chromosomes.contains_key(gene.0) {
                self.insert_chromosome(&gene.0, Chromosome::new());
            }

            self.chromosomes.get_mut(gene.0).unwrap().insert_gene(&gene.1, gene.2);
        }
    }

    /// Sets the value of the specified chromosome.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let mut genome = Genome::new(vec![("examples", "example1", Gene::default())]);
    /// let _ = genome.set_chromosome("examples", Chromosome::default());
    /// assert!(genome.gene("examples", "example1").is_none());
    /// ```
    #[inline]
    pub fn set_chromosome(&mut self, chromo_name: &str, chromo: Chromosome) -> Result<(), &str> {
        *self.chromosomes.get_mut(chromo_name).ok_or("No chromosome by that name exists")? = chromo;

        Ok(())
    }

    /// Sets the value of the specified gene.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use genetic_optimization::prelude::*;
    /// 
    /// let mut genome = Genome::new(vec![("examples", "example1", Gene::new(1.0))]);
    /// let _ = genome.set_gene("examples", "example1", Gene::default());
    /// assert_eq!(genome.gene("examples", "example1").unwrap().value(), 0.0);
    /// ```
    #[inline]
    pub fn set_gene(&mut self, chromo_name: &str, gene_name: &str, gene: Gene) -> Result<(), &str> {
        self.chromosomes.get_mut(chromo_name).ok_or("No chromosome by that name exists")?.set_gene(gene_name, gene)
    }

    /// Saves `self` to file.
    #[inline]
    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), &str> {
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
        let mut display_obj: HashMap<String, HashMap<String, f32>> = HashMap::new();
        for chromo in &self.chromosomes {
            display_obj.insert(chromo.0.to_string(), HashMap::new());
            for gene in chromo.1.genes() {
                display_obj.get_mut(chromo.0).unwrap().insert(gene.0.to_string(), gene.1.value());
            }
        }
        write!(f, "{:#?}", display_obj)
    }
}