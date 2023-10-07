type Float = f32;
const G: Float = 1.184e-4;
const ROCK_DENSITY: Float = 2.803e12; //5 g/cm^3 = 5000 kg/m^3
const ATM_DENSITY: Float = 5.61e10; //0.1 g/cm^3, a rough estimate
const DIMENSIONALITY: usize = 2;
const ROCK_MASS_TO_VOLUME_FACTOR: Float = 3. / (4. * std::f32::consts::PI * ROCK_DENSITY);
const ATM_MASS_TO_VOLUME_FACTOR: Float = 3. / (4. * std::f32::consts::PI * ATM_DENSITY);

pub(crate) struct Body {
    position: Vec<Float>,
    velocity: Vec<Float>,
    rocky_mass: Float,
    gas_mass: Float,
}

impl Body {
    fn new(position: Vec<Float>, velocity: Vec<Float>, rocky_mass: Float, gas_mass: Float) -> Body {
        Body {
            position,
            velocity,
            rocky_mass,
            gas_mass,
        }
    }

    fn rocky_radius(&self) -> Float {
        (self.rocky_mass * ROCK_MASS_TO_VOLUME_FACTOR).cbrt()
    }

    fn atm_radius(&self) -> Float {
        let rocky_radius = self.rocky_radius();
        (self.gas_mass * ATM_MASS_TO_VOLUME_FACTOR + rocky_radius.powi(3)).cbrt() - rocky_radius
    }

    fn radius(&self) -> Float {
        self.rocky_radius() + self.atm_radius()
    }
}
