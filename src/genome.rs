use std::collections::HashMap;

use crate::prelude::*;

#[derive(Clone, Debug, Default, PartialEq)]
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
    /// let genome = Genome::new(vec!["chromosome_label", "gene_label", Gene::default()]);
    /// ```
    #[inline]
    pub fn new(genes: Vec<(&str, &str, Gene)>) -> Self {
        let mut genome = Self { 
            ..Default::default()
        };

        genome.insert_genes(genes);

        genome
    }

    #[inline]
    /// Returns a reference to `self.chromosomes`.
    pub fn chromosomes(&self) -> &HashMap<String, Chromosome> {
        &self.chromosomes
    }

    #[inline]
    /// Returns a list of all genes, their labels, and their containers' labels.
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
    /// let mut genome = Genome::new(vec!["examples", "example1", Gene::default()]);
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
    /// let mut genome = Genome::new(vec!["examples", "example1", Gene::default()]);
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
    /// genome.insert_genes(vec!["chromosome_label", "gene_label", the_gene_itself]);
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

    /// Sets the value of the given gene.
    #[inline]
    pub fn set_gene(&mut self, chromo_name: &str, gene_name: &str, gene: Gene) {
        self.chromosomes.get_mut(chromo_name).unwrap().set_gene(gene_name, gene);
    }
}