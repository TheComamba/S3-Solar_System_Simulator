use crate::initial_parameters::{Float, DIMENSIONALITY, ROCK_DENSITY};
use rand_distr::{Distribution, Normal};

pub(crate) struct Body {
    position: Vec<Float>,
    velocity: Vec<Float>,
    mass: Float,
}

impl Body {
    fn random_vector(variance: Float) -> Vec<Float> {
        let mut vector = vec![0.; DIMENSIONALITY];
        let distribution = Normal::new(0., variance).unwrap();
        for i in 0..DIMENSIONALITY {
            vector[i] = distribution.sample(&mut rand::thread_rng());
        }
        vector
    }

    pub(crate) fn new(position_variance: Float, velocity_variance: Float, mass: Float) -> Body {
        Body {
            position: Self::random_vector(position_variance),
            velocity: Self::random_vector(velocity_variance),
            mass,
        }
    }

    pub(crate) fn radius(&self) -> Float {
        const MASS_TO_VOLUME_FACTOR: Float = 3. / (4. * std::f32::consts::PI * ROCK_DENSITY);
        (self.mass * MASS_TO_VOLUME_FACTOR).cbrt()
    }
}
