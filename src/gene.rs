use rand::Rng;

/// A specific parameter that will be optimized by the simulation.
#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct Gene {
    value: f32,
    min: f32,
    max: f32,
}

impl Gene {
    /// Creates a new `Gene`.
    #[inline]
    pub fn new(value: f32) -> Self {
        Self { 
            value,
            min: 0.0,
            max: 1.0,
        }
    }

    /// Creates a new `Gene` with the specified min and max values.
    #[inline]
    pub fn new_with_range(value: f32, min: f32, max: f32) -> Self {
        Self { 
            value,
            min,
            max,
        }
    }

    /// Returns `self.value`.
    #[inline]
    pub fn value(&self) -> f32 {
        self.value
    }

    /// Returns `self.min`.
    #[inline]
    pub fn min(&self) -> f32 {
        self.min
    }

    /// Returns `self.max`.
    #[inline]
    pub fn max(&self) -> f32 {
        self.max
    }

    /// Creates a mutated `Gene` from `self`.
    #[inline]
    pub fn mutate(&self, randomness_weight: f32) -> Gene {
        let mut rng = rand::thread_rng();
        let mut gene = self.clone();

        gene.value = 1.0 / (randomness_weight + 1.0) * (self.value + randomness_weight * rng.gen_range(self.min..self.max));

        gene
    }
}