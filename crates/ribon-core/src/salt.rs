//! Monovalent-salt corrections for the Turner nearest-neighbor model.
//!
//! The implementation follows the published polyelectrolyte expressions and
//! returns corrections in centi-kcal/mol.

use std::f64::consts::PI;

pub(crate) const STANDARD_MOLAR: f64 = 1.021;
const GAS_CONSTANT_CAL: f64 = 1.987_17;
const RODS_DISTANCE: f64 = 20.0;
const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;
const BACKBONE_LENGTH: f64 = 6.0;
const HELICAL_RISE: f64 = 2.8;

#[derive(Clone, Copy, Debug)]
pub(crate) struct SaltCorrections {
    pub(crate) molar: f64,
    pub(crate) stack_centi: i32,
    pub(crate) ml_base_centi: i32,
    pub(crate) ml_closing_centi: i32,
    pub(crate) duplex_init_centi: i32,
}

impl SaltCorrections {
    pub(crate) fn new(molar: f64, temperature_celsius: f64) -> Self {
        if (molar - STANDARD_MOLAR).abs() < 1.0e-12 {
            return Self {
                molar,
                stack_centi: 0,
                ml_base_centi: 0,
                ml_closing_centi: 0,
                duplex_init_centi: 0,
            };
        }
        let temperature_kelvin = temperature_celsius + 273.15;
        let stack_centi = round_centi(stack(molar, temperature_kelvin, HELICAL_RISE));
        let values: Vec<_> = (0..=24)
            .map(|length| loop_raw(length, molar, temperature_kelvin, BACKBONE_LENGTH))
            .collect();
        let (ml_base_centi, ml_closing_centi) = linear_fit(&values, 6, 24);
        let duplex_init_centi = round_centi((molar / STANDARD_MOLAR).ln() * -45.324);
        Self {
            molar,
            stack_centi,
            ml_base_centi,
            ml_closing_centi,
            duplex_init_centi,
        }
    }

    pub(crate) fn loop_centi(&self, length: usize, temperature_celsius: f64) -> i32 {
        if (self.molar - STANDARD_MOLAR).abs() < 1.0e-12 || length == 0 {
            0
        } else {
            round_centi(loop_raw(
                length,
                self.molar,
                temperature_celsius + 273.15,
                BACKBONE_LENGTH,
            ))
        }
    }
}

fn epsilon_r(temperature: f64) -> f64 {
    5321.0 / temperature + 233.76 - 0.9297 * temperature
        + 1.417 * temperature * temperature / 1000.0
        - 0.8292 * temperature * temperature * temperature / 1_000_000.0
}

fn bjerrum_length(temperature: f64) -> f64 {
    167_100.052 / (temperature * epsilon_r(temperature))
}

fn kappa(molar: f64, temperature: f64) -> f64 {
    (bjerrum_length(temperature) * molar).sqrt() / 8.1284
}

fn tau_ds(temperature: f64, helical_rise: f64) -> f64 {
    (1.0 / helical_rise).min(1.0 / bjerrum_length(temperature))
}

fn tau_ss(temperature: f64, backbone_length: f64) -> f64 {
    (1.0 / backbone_length).min(1.0 / bjerrum_length(temperature))
}

fn pairing_constant(temperature: f64, helical_rise: f64) -> f64 {
    let tau = tau_ds(temperature, helical_rise);
    2.0 * (GAS_CONSTANT_CAL / 1000.0)
        * temperature
        * bjerrum_length(temperature)
        * helical_rise
        * tau
        * tau
}

fn screening_hyper_term(y: f64) -> f64 {
    let a = 1.0 / (y.powi(6) / (2.0 * PI).powi(6) + 1.0);
    let b = y.powi(4) / (36.0 * PI.powi(4)) - y.powi(3) / (24.0 * PI * PI)
        + y * y / (2.0 * PI * PI)
        - y / 2.0;
    let c = (2.0 * PI / y).ln() - 1.963_51;
    a * b + (1.0 - a) * c
}

fn loop_aux(kappa_length: f64, length: usize, temperature: f64, backbone: f64) -> f64 {
    let tau = tau_ss(temperature, backbone);
    let a = (GAS_CONSTANT_CAL / 1000.0)
        * temperature
        * bjerrum_length(temperature)
        * length as f64
        * backbone
        * tau
        * tau;
    let b = kappa_length.ln() - (PI / 2.0).ln()
        + EULER_MASCHERONI
        + screening_hyper_term(kappa_length)
        + (1.0 - (-kappa_length).exp() + kappa_length * exponential_integral_e1(kappa_length))
            / kappa_length;
    a * b * 100.0
}

fn loop_raw(length: usize, molar: f64, temperature: f64, backbone: f64) -> f64 {
    if length == 0 {
        return 0.0;
    }
    let span = length as f64 * backbone;
    loop_aux(
        kappa(molar, temperature) * span,
        length,
        temperature,
        backbone,
    ) - loop_aux(
        kappa(STANDARD_MOLAR, temperature) * span,
        length,
        temperature,
        backbone,
    )
}

fn stack(molar: f64, temperature: f64, helical_rise: f64) -> f64 {
    100.0
        * pairing_constant(temperature, helical_rise)
        * (bessel_k0(RODS_DISTANCE * kappa(molar, temperature))
            - bessel_k0(RODS_DISTANCE * kappa(STANDARD_MOLAR, temperature)))
}

fn linear_fit(values: &[f64], lower: usize, upper: usize) -> (i32, i32) {
    let count = (upper - lower + 1) as f64;
    let mut sum_x = 0.0;
    let mut sum_xx = 0.0;
    let mut sum_y = 0.0;
    let mut sum_xy = 0.0;
    for (index, &value) in values.iter().enumerate().take(upper + 1).skip(lower) {
        let x = index as f64;
        sum_x += x;
        sum_xx += x * x;
        sum_y += value;
        sum_xy += x * value;
    }
    let denominator = count * sum_xx - sum_x * sum_x;
    let slope = (count * sum_xy - sum_x * sum_y) / denominator;
    let intercept = (sum_y * sum_xx - sum_x * sum_xy) / denominator;
    (round_centi(slope), round_centi(intercept))
}

fn round_centi(value: f64) -> i32 {
    if value < 0.0 {
        (value - 0.5) as i32
    } else {
        (value + 0.5) as i32
    }
}

// Convergent continued fraction / power series for E1(x), x > 0. Both paths
// stop from f64 convergence or representational stagnation, never from a
// hidden iteration ceiling.
fn exponential_integral_e1(x: f64) -> f64 {
    const EULER: f64 = 0.577_215_664_901_532_9;
    const EPSILON: f64 = 1.0e-15;
    if x > 1.0 {
        let mut b = x + 1.0;
        let mut c = 1.0 / f64::MIN_POSITIVE;
        let mut d = 1.0 / b;
        let mut h = d;
        let mut index = 1usize;
        loop {
            let a = -(index as f64).powi(2);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            let delta = c * d;
            let previous = h;
            h *= delta;
            if (delta - 1.0).abs() < EPSILON || h == previous {
                break;
            }
            index += 1;
        }
        h * (-x).exp()
    } else {
        let mut answer = -x.ln() - EULER;
        let mut factorial = 1.0;
        let mut index = 1usize;
        loop {
            factorial *= -x / index as f64;
            let delta = -factorial / index as f64;
            let previous = answer;
            answer += delta;
            if delta.abs() < answer.abs().max(1.0) * EPSILON || answer == previous {
                break;
            }
            index += 1;
        }
        answer
    }
}

// Directly convergent power series for I0. The salt model only needs I0 in
// the K0 small-argument expansion (x <= 2), where this series is rapid.
fn bessel_i0(x: f64) -> f64 {
    let y = x * x / 4.0;
    let mut term = 1.0;
    let mut sum = 1.0;
    let mut index = 1usize;
    loop {
        term *= y / (index * index) as f64;
        let previous = sum;
        sum += term;
        if term.abs() <= sum.abs() * f64::EPSILON || sum == previous {
            return sum;
        }
        index += 1;
    }
}

fn bessel_k0(x: f64) -> f64 {
    if x <= 2.0 {
        let y = x * x / 4.0;
        let mut term = 1.0;
        let mut harmonic = 0.0;
        let mut series = 0.0;
        let mut index = 1usize;
        loop {
            harmonic += 1.0 / index as f64;
            term *= y / (index * index) as f64;
            let delta = harmonic * term;
            let previous = series;
            series += delta;
            if delta.abs() <= series.abs().max(1.0) * f64::EPSILON || series == previous {
                break;
            }
            index += 1;
        }
        -((x / 2.0).ln() + EULER_MASCHERONI) * bessel_i0(x) + series
    } else {
        bessel_k0_integral(x)
    }
}

/// K0(x) = integral_0^infinity exp(-x cosh(t)) dt, evaluated by adaptive
/// Simpson quadrature. The upper tail is extended until its convex-exponential
/// bound is below the requested f64-relative tolerance; subdivision ends only
/// on its error estimate or when no further midpoint is representable.
fn bessel_k0_integral(x: f64) -> f64 {
    let integrand = |t: f64| (-x * t.cosh()).exp();
    let relative_tolerance = 32.0 * f64::EPSILON;
    let mut upper = 1.0f64;
    loop {
        let tail_bound = integrand(upper) / (x * upper.sinh());
        if tail_bound <= relative_tolerance * integrand(0.0).max(f64::MIN_POSITIVE)
            || tail_bound == 0.0
        {
            break;
        }
        upper *= 2.0;
    }
    let fa = integrand(0.0);
    let fm = integrand(upper / 2.0);
    let fb = integrand(upper);
    let whole = upper * (fa + 4.0 * fm + fb) / 6.0;
    let tolerance = relative_tolerance * whole.abs().max(f64::MIN_POSITIVE);
    let mut stack = vec![(0.0, upper, fa, fm, fb, whole, tolerance)];
    let mut result = 0.0;
    while let Some((left, right, f_left, f_middle, f_right, parent, tolerance)) = stack.pop() {
        let middle = (left + right) / 2.0;
        let left_middle = (left + middle) / 2.0;
        let right_middle = (middle + right) / 2.0;
        let f_left_middle = integrand(left_middle);
        let f_right_middle = integrand(right_middle);
        let left_integral = (middle - left) * (f_left + 4.0 * f_left_middle + f_middle) / 6.0;
        let right_integral = (right - middle) * (f_middle + 4.0 * f_right_middle + f_right) / 6.0;
        let refined = left_integral + right_integral;
        let error = (refined - parent).abs();
        if error <= 15.0 * tolerance
            || left_middle == left
            || left_middle == middle
            || right_middle == middle
            || right_middle == right
        {
            result += refined + (refined - parent) / 15.0;
        } else {
            stack.push((
                middle,
                right,
                f_middle,
                f_right_middle,
                f_right,
                right_integral,
                tolerance / 2.0,
            ));
            stack.push((
                left,
                middle,
                f_left,
                f_left_middle,
                f_middle,
                left_integral,
                tolerance / 2.0,
            ));
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_salt_is_exactly_neutral() {
        let correction = SaltCorrections::new(STANDARD_MOLAR, 37.0);
        assert_eq!(correction.stack_centi, 0);
        assert_eq!(correction.ml_base_centi, 0);
        assert_eq!(correction.ml_closing_centi, 0);
        assert_eq!(correction.duplex_init_centi, 0);
        assert_eq!(correction.loop_centi(12, 37.0), 0);
    }

    #[test]
    fn special_functions_are_finite_over_model_range() {
        for value in [0.01, 0.1, 1.0, 2.0, 10.0, 100.0] {
            assert!(bessel_k0(value).is_finite());
            assert!(exponential_integral_e1(value).is_finite());
            assert!(exponential_integral_e1(value) > 0.0);
        }
    }

    #[test]
    fn k0_matches_high_precision_reference_values() {
        for (argument, expected) in [
            (0.1, 2.427_069_024_702_016_4),
            (1.0, 0.421_024_438_240_708_3),
            (2.0, 0.113_893_872_749_533_44),
            (10.0, 1.778_006_231_616_765e-5),
            (100.0, 4.656_628_229_175_902e-45),
        ] {
            let actual = bessel_k0(argument);
            assert!(
                ((actual - expected) / expected).abs() < 2.0e-13,
                "K0({argument}) = {actual:.17e}, expected {expected:.17e}"
            );
        }
    }
}
