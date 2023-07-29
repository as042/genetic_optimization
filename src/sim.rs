use rand::Rng;

use crate::prelude::*;

impl Genome {
    #[inline]
    pub fn simulate(&self, iterations: usize, evaluator: fn(Genome) -> f32) -> Genome {
        let init_species = gen_random_species(self, 16, 1.0);

        let mut best_genome = self.clone();
        for k in 0..iterations {
            
        }
        
        best_genome
    }
}

#[inline]
fn gen_random_species(default: &Genome, num: usize, range: f32) -> Vec<Genome> {
    let mut rng = rand::thread_rng();

    let mut species = Vec::default();
    for _ in 0..num {
        let mut spec = default.clone();
        for gene in default.genes() {
            spec.set_gene(gene.0, gene.1, Gene::new(gene.2.value + range * rng.gen::<f32>()));
        }

        species.push(spec)
    }

    species
}