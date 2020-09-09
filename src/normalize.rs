use std::f32;
use std::vec::Vec;
const LOWER: f32 = -3.;
const UPPER: f32 = 3.;
const DX: f32 = 0.001;
/// Enum for normalization
#[derive(Debug, Clone, Copy)]
pub enum NormalizeType {
    /// Z-normalization
    Z,
    ///Max min normalization
    MaxMin,
}
fn z_normalize(xs: &[f32]) -> Vec<f32> {
    let len = xs.len() as f64;
    let (sum, sqsum): (f64, f64) = xs.iter().fold((0., 0.), |(sum, sqsum), &x| {
        let x = x as f64;
        (sum + x, sqsum + x * x)
    });
    let mean = sum / len;
    let stdev = (sqsum / len - mean * mean).sqrt();
    // if stdev <= 0.000001 {
    //     eprintln!("varianve is strictly zero!");
    // }
    let res = xs
        .iter()
        .map(|&x| ((x as f64 - mean) / stdev) as f32)
        .collect::<Vec<_>>();
    res
}
#[test]
fn test() {
    let v = (0..100).map(|e| f32::from((e - 50 as i8))).collect();
    z_normalize(&v);
    assert!(true);
}
fn max_min_normalize(xs: &[f32]) -> Vec<f32> {
    let (max, min) = xs.iter().fold((f32::MIN, f32::MAX), |(big, small), &x| {
        (big.max(x), small.min(x))
    });
    if (max - min).abs() < std::f32::EPSILON {
        panic!("the input is allwhere the same !");
    }
    xs.iter().map(|x| (x - min) / (max - min)).collect()
}
/// Normalize vector by given normalized type
pub fn normalize(xs: &[f32], mode: NormalizeType) -> Vec<f32> {
    match mode {
        NormalizeType::MaxMin => max_min_normalize(xs),
        NormalizeType::Z => z_normalize(xs),
    }
}

fn z_normalize_mut(xs: &mut [f32]) {
    let len = xs.len() as f32;
    let (sum, sqsum): (f32, f32) = xs
        .iter()
        .fold((0., 0.), |(sum, sqsum), &x| (sum + x, sqsum + x * x));
    let mean = sum / len;
    let stdev = (sqsum / len - mean * mean).sqrt();
    // if stdev <= 0.000001 {
    //     eprintln!("varianve is strictly zero!");
    // }
    for x in xs.iter_mut() {
        *x = (*x - mean) / stdev;
    }
}
fn max_min_normalize_mut(xs: &mut [f32]) {
    let max = xs
        .iter()
        .fold(-1000000., |max, &x| if max < x { x } else { max });
    let min = xs
        .iter()
        .fold(1000000., |min, &x| if min < x { min } else { x });
    for x in xs.iter_mut() {
        *x = (*x - min) / max;
    }
}

/// Normalize the given mutable refernece.
pub fn normalize_mut(xs: &mut [f32], mode: NormalizeType) {
    match mode {
        NormalizeType::Z => z_normalize_mut(xs),
        NormalizeType::MaxMin => max_min_normalize_mut(xs),
    }
}

/// histgram equization.
/// This function assumes that the argument is in the range of [-3,3].
/// If not, return error.
/// If the argument satisfy the condition, this function compute the equilized
/// vector, and returns the result and its conversion function.
/// When you need to apply the same histgram conversion to other vector,
/// use histgram_modify_by(), but make sure that is not for "equalizing".
pub fn histgram_equalization(xs: &Vec<f32>) -> (Vec<f32>, Vec<f32>) {
    let cdf = cumulative_dist_function(xs);
    (histgram_modify(xs, &cdf), cdf)
}

fn cumulative_dist_function(xs: &Vec<f32>) -> Vec<f32> {
    let n = xs.len() as f32;
    let mut xs = xs.clone();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut xs = xs.iter().peekable();
    let size = ((UPPER - LOWER) / DX) as usize;
    let mut cdf = Vec::with_capacity(size);
    let mut current_bound = LOWER;
    let mut acm_freq = 0.;
    for _ in 0..size {
        // calculate the number of element which is smaller than current bound
        loop {
            if let Some(&&x) = xs.peek() {
                if x < current_bound {
                    acm_freq += 1.;
                    xs.next();
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        current_bound += DX;
        cdf.push(acm_freq / n);
    }
    acm_freq += xs.count() as f32;
    cdf[size - 1] = acm_freq / n;
    cdf
}

/// histgram equilization by using given cumulative distribution function
pub fn histgram_modify(xs: &[f32], cdf: &[f32]) -> Vec<f32> {
    let n = xs.len();
    let mut res = vec![0.; n];
    let size = cdf.len();
    let range = UPPER - LOWER;
    for k in 0..n {
        let i = ((xs[k] - LOWER) / DX).floor() as usize;
        if i >= size - 1 {
            res[k] = cdf[size - 1] * range + LOWER;
        } else if i == 0 {
            res[k] = cdf[0] * range + LOWER;
        } else {
            res[k] = (cdf[i + 1] - cdf[i]) * range / DX * (xs[k] - LOWER - i as f32 * DX)
                + cdf[i] * range
                + LOWER;
        }
    }
    res
}
