use std::collections::HashMap;
use serde::{Serialize, Deserialize};

use crate::prelude::*;

/// A list of semi-related genes.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Chromosome {
    genes: HashMap<String, Gene>,
}

impl Chromosome {
    /// Creates a new `Chromosome`.
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a gene to `self` and returns `self`.
    #[inline]
    pub fn add_gene(&mut self, gene_name: impl Into<String>, gene: Gene) -> &mut Self {
        self.insert_gene(gene_name, gene);
        self
    }

    /// Builds `self` and returns an owned value.
    #[inline]
    pub fn build(&mut self) -> Self {
        self.to_owned()
    }

    /// Returns the specified gene, if it exists.
    #[inline]
    pub fn gene(&self, gene_name: impl Into<String>) -> Option<&Gene> {
        self.genes.get(&gene_name.into())
    }

    /// Returns a reference to `self.genes`.
    #[inline]
    pub fn genes(&self) -> &HashMap<String, Gene> {
        &self.genes
    }

    /// Inserts a `Gene` into `self`.
    #[inline]
    pub fn insert_gene(&mut self, gene_name: impl Into<String>, gene: Gene) {
        self.genes.insert(gene_name.into(), gene);
    }

    /// Inserts all given genes into `self`.
    #[inline]
    pub fn insert_genes(&mut self, genes: Vec<(impl Into<String>, Gene)>) {
        for gene in genes {
            self.insert_gene(gene.0, gene.1);
        }
    }

    /// Sets the value of the specified gene.
    #[inline]
    pub fn set_gene(&mut self, gene_name: impl Into<String>, gene: Gene) -> Result<(), String> {
        *self.genes.get_mut(&gene_name.into()).ok_or("No gene by that name exists")? = gene;

        Ok(())
    }
}