use std::collections::HashMap;

use crate::prelude::*;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Chromosome {
    genes: HashMap<String, Gene>,
}

impl Chromosome {
    /// Creates a new `Chromosome`.
    #[inline]
    pub fn new() -> Self {
        Self { 
            ..Default::default()
        }
    }

    #[inline]
    pub fn genes(&self) -> &HashMap<String, Gene> {
        &self.genes
    }

    /// Inserts a `Gene` into `self`.
    #[inline]
    pub fn insert_gene(&mut self, gene_name: &str, gene: Gene) {
        self.genes.insert(gene_name.to_string(), gene);
    }

    /// Inserts all given genes into `self`.
    #[inline]
    pub fn insert_genes(&mut self, genes: Vec<(&str, Gene)>) {
        for gene in genes {
            self.insert_gene(gene.0, gene.1);
        }
    }

    /// Sets the value of the given gene.
    #[inline]
    pub fn set_gene(&mut self, gene_name: &str, gene: Gene) {
        *self.genes.get_mut(gene_name).unwrap() = gene;
    }
}