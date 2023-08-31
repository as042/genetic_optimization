use rand::Rng;
use serde::{Serialize, Deserialize};

/// A specific parameter that will be optimized by the simulation.
#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Gene {
    value: f64,
    min: f64,
    max: f64,
}

impl Gene {
    /// Creates a new `Gene`.
    #[inline]
    pub fn new(value: f64) -> Self {
        Self { 
            value,
            min: 0.0,
            max: 1.0,
        }
    }

    /// Creates a new `Gene` with the specified min and max values.
    #[inline]
    pub fn new_with_range(value: f64, min: f64, max: f64) -> Self {
        Self { 
            value,
            min,
            max,
        }
    }

    /// Returns `self.value`.
    #[inline]
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Returns `self.min`.
    #[inline]
    pub fn min(&self) -> f64 {
        self.min
    }

    /// Returns `self.max`.
    #[inline]
    pub fn max(&self) -> f64 {
        self.max
    }

    /// Creates a mutated `Gene` from `self`.
    #[inline]
    pub fn mutate(&self, randomness_weight: f64) -> Gene {
        let mut rng = rand::thread_rng();
        let mut gene = self.clone();

        gene.value = 1.0 / (randomness_weight + 1.0) * (self.value + randomness_weight * rng.gen_range(self.min..self.max));

        gene
    }
}