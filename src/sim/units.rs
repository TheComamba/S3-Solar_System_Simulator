pub(crate) type Float = f32;
pub(crate) const DIMENSIONALITY: usize = 2;
//pub(crate) const PI: Float = std::f64::consts::PI as Float;
pub(crate) const PI: Float = 3.14159265358979323846264338327950288;
//Defined:
//Measuring masses in earth masses
pub(crate) const MASS_TO_KG: Float = 5.972e24;
//Measuring distances in AU
pub(crate) const DISTANCE_TO_M: Float = 1.496e11;
//Measuring time in years
pub(crate) const TIME_TO_SECONDS: Float = 3.154e7;

//Derived:
// AU = (G M year^2 / (4 Pi^2))^(1/3), where M = Mass of sun = 333000 M_E
// => G = 4 Pi^2 AU^3 / (333000 M_E year^2)
pub(crate) const G: Float = 4. * PI * PI / 333_000.;
// Assuming 5 g/cm^3 = 5000 kg/m^3
pub(crate) const ROCK_DENSITY: Float =
    5_000. * DISTANCE_TO_M * DISTANCE_TO_M * DISTANCE_TO_M / MASS_TO_KG; //
pub(crate) const VELOCITY_TO_KM_PER_S: Float = DISTANCE_TO_M / TIME_TO_SECONDS / 1e3;
