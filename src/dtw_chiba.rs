use num::Float;
#[inline]
fn dp_to_matrix_idx(x: usize, band: usize) -> (usize, usize) {
    //    println!("{}->({},{})",x,x/band,x/band + x % band - (band-1)/2);
    (x / band, x / band + x % band - (band - 1) / 2)
}
#[test]
fn dp_to_matrix_test() {
    assert_eq!(dp_to_matrix_idx(4, 3), (1, 1));
    assert_eq!(dp_to_matrix_idx(10, 3), (3, 3));
    assert_eq!(dp_to_matrix_idx(22, 3), (7, 7));
    assert_eq!(dp_to_matrix_idx(22, 5), (4, 4));
    assert_eq!(dp_to_matrix_idx(26, 5), (5, 4));
    assert_eq!(dp_to_matrix_idx(4, 5), (0, 2));
}
#[inline]
fn matrix_to_dp_idx(i: usize, j: usize, band: usize) -> usize {
    //    println!("({},{}) -> {}",i,j,band*i + j - i + band/2);
    band * i + j - i + band / 2
}
#[test]
fn matrix_to_dp_test() {
    assert_eq!(matrix_to_dp_idx(0, 0, 3), 1);
    assert_eq!(matrix_to_dp_idx(7, 7, 3), 22);
    assert_eq!(matrix_to_dp_idx(4, 3, 3), 12);
    assert_eq!(matrix_to_dp_idx(5, 5, 5), 27);
    assert_eq!(matrix_to_dp_idx(3, 5, 5), 19);
    assert_eq!(matrix_to_dp_idx(1, 3, 5), 9);
}

#[test]
fn inversivity() {
    let band = 11;
    let n = 20;
    for i in 0..n {
        for j in 0..n {
            if j < i + band / 2 && i < j + band / 2 {
                let dp = matrix_to_dp_idx(i, j, band);
                let (n, m) = dp_to_matrix_idx(dp, band);
                assert_eq!((i, j), (n, m));
            }
        }
    }
}

#[test]
fn inversivity2() {
    let band = 21;
    let n = 30;
    for dp in 0..band * n {
        if dp_valid_idx(dp, band, n) {
            let (i, j) = dp_to_matrix_idx(dp, band);
            let m = matrix_to_dp_idx(i, j, band);
            debug_assert!(dp == m, "{} -> ({},{}) -> {}", dp, i, j, m);
        }
    }
}
#[inline]
fn dp_valid_idx(x: usize, band: usize, len: usize) -> bool {
    x / band + x % band <= (band - 1) / 2 + len && x / band + x % band > (band - 1) / 2
}

fn nextpos<D, F, T>(
    x1: &[D],
    x2: &[D],
    dp: &[T],
    x1pos: usize,
    x2pos: usize,
    dist: &F,
    band: usize,
) -> (usize, usize)
where
    F: Fn(&D, &D) -> T,
    T: Float,
{
    let current_score =
        dp[matrix_to_dp_idx(x1pos, x2pos, band)] - dist(&x1[x1pos - 1], &x2[x2pos - 1]);
    let match_path = dp[matrix_to_dp_idx(x1pos - 1, x2pos - 1, band)] - current_score;
    let gap_path = dp[matrix_to_dp_idx(x1pos, x2pos - 1, band)] - current_score;
    let del_path = dp[matrix_to_dp_idx(x1pos - 1, x2pos, band)] - current_score;
    if match_path < gap_path && match_path < del_path {
        (x1pos - 1, x2pos - 1)
    } else if gap_path < match_path && gap_path < del_path {
        (x1pos, x2pos - 1)
    } else {
        (x1pos - 1, x2pos)
    }
}

#[inline]
fn detect_err_about_band(x1len: usize, x2len: usize, band: usize) -> Result<(), String> {
    match band {
        n if n == 1 => Err("Band width should be greater than or equal to 3".to_string()),
        n if n % 2 == 0 => Err("Band width should be an odd number".to_string()),
        n if n > x1len || n > x2len => {
            Err("Band width should be smaller than query length".to_string())
        }
        _ if x1len != x2len => Err("The two sequence is not the same length".to_string()),
        _ => Ok(()),
    }
}

#[test]
fn err_detect1() {
    let x1len = 10;
    let x2len = 10;
    let band = 3;
    assert!(detect_err_about_band(x1len, x2len, band).is_ok());
}
#[test]
fn err_detect2() {
    let x1len = 10;
    let x2len = 10;
    let band = 1;
    assert!(detect_err_about_band(x1len, x2len, band).is_err());
}
#[test]
fn err_detect3() {
    let x1len = 10;
    let x2len = 10;
    let band = 4;
    assert!(detect_err_about_band(x1len, x2len, band).is_err());
}
#[test]
fn err_detect4() {
    let x1len = 10;
    let x2len = 10;
    let band = 11;
    assert!(detect_err_about_band(x1len, x2len, band).is_err());
}
#[test]
fn err_detect5() {
    let x1len = 201;
    let x2len = 200;
    let band = 100;
    assert!(detect_err_about_band(x1len, x2len, band).is_err());
}

#[allow(dead_code)]
fn dpflush<T>(dp: &Vec<T>, band: usize, n: usize)
where
    T: Float,
{
    println!();
    for i in 0..n + 1 {
        for j in 0..n + 1 {
            if j <= i + band / 2 && i <= j + band / 2 {
                let x = dp[matrix_to_dp_idx(i, j, band)];
                if x.is_infinite() {
                    print!("iii");
                } else {
                    print!("{:3.0}", dp[matrix_to_dp_idx(i, j, band)].to_f32().unwrap());
                };
            } else {
                print!("eee");
            }
        }
        println!();
    }
    // for i in 0..dp.len(){
    //     if i % band == 0{
    //         println!();
    //     }
    //     print!("{:3}",dp[i].to_f32().unwrap());
    // }
}

type TRACE = (Vec<usize>, Vec<usize>);
pub fn dtw_chiba<D, F, T>(
    x1: &[D],
    x2: &[D],
    dist: &F,
    band: usize,
) -> Result<(T, TRACE, usize), String>
where
    F: Fn(&D, &D) -> T,
    T: Float,
{
    detect_err_about_band(x1.len(), x2.len(), band)?;
    let n = x1.len();
    let mut dp = vec![num::zero(); (n + 1) * band];
    for i in 0..(band + 1) / 2 {
        dp[matrix_to_dp_idx(i, 0, band)] = Float::infinity();
    }
    for j in 0..(band + 1) / 2 {
        dp[matrix_to_dp_idx(0, j, band)] = Float::infinity();
    }
    for i in (band + 1) / 2..n + 1 {
        dp[matrix_to_dp_idx(i, i - (band - 1) / 2, band)] = Float::infinity();
    }
    for j in (band + 1) / 2..n + 1 {
        dp[matrix_to_dp_idx(j - (band - 1) / 2, j, band)] = Float::infinity();
    }
    dp[matrix_to_dp_idx(0, 0, band)] = num::zero();
    // dynamic time warping
    //dpflush(&dp,band,n);
    for idx in (band / 2) + band..(n + 1) * band {
        let current: T = dp[idx];
        if !current.is_infinite() && dp_valid_idx(idx, band, n) {
            let (i, j) = dp_to_matrix_idx(idx, band);
            let prev1 = matrix_to_dp_idx(i - 1, j - 1, band);
            let prev2 = matrix_to_dp_idx(i, j - 1, band);
            let prev3 = matrix_to_dp_idx(i - 1, j, band);
            let cost = dist(&x1[i - 1], &x2[j - 1]);
            let min: T = if dp[prev2] < dp[prev3] && dp[prev2] < dp[prev1] {
                dp[prev2]
            } else if dp[prev3] < dp[prev2] && dp[prev3] < dp[prev1] {
                dp[prev3]
            } else {
                dp[prev1]
            };
            dp[idx] = if !min.is_infinite() { cost + min } else { min };
        }
    }
    //dpflush(&dp,band,n);
    // determine warping path
    let mut x1path = vec![];
    let mut x2path = vec![];
    let mut x1pos = n;
    let mut x2pos = n;
    while x1pos > 0 {
        x1path.push(x1pos - 1);
        x2path.push(x2pos - 1);
        let next_pos = nextpos(x1, x2, &dp, x1pos, x2pos, dist, band);
        x1pos = next_pos.0;
        x2pos = next_pos.1;
    }
    x1path.reverse();
    x2path.reverse();
    Ok((dp[n * band + band / 2], (x1path, x2path), 0))
}

#[test]
fn phony_test() {
    let x1 = vec![1, 2, 3];
    let x2 = vec![1, 2, 3];
    let d = |x: &i32, y: &i32| (x - y).abs() as f32;
    let band = 3;
    let (score, (xpath, ypath), start) = dtw_chiba(&x1, &x2, &d, band).unwrap();
    assert_eq!(score, 0.0);
    assert_eq!(xpath, vec![0, 1, 2]);
    assert_eq!(ypath, vec![0, 1, 2]);
    assert_eq!(start, 0);
}

#[test]
fn too_large_band() {
    let mut x1: Vec<_> = (0..10).collect();
    x1.push(0);
    x1.push(0);
    let mut x2: Vec<_> = vec![0];
    x2.append(&mut (0..10).collect());
    x2.push(0);
    let d = |x: &i32, y: &i32| (x - y).abs() as f32;
    let band = 21;
    assert!(dtw_chiba(&x1, &x2, &d, band).is_err());
}

#[test]
fn flat_test() {
    let x1 = vec![0; 20];
    let x2 = vec![0; 20];
    let d = |x: &i32, y: &i32| (x - y).abs() as f32;
    let band = 11;
    let (score, _, _) = dtw_chiba(&x1, &x2, &d, band).unwrap();
    assert_eq!(score, 0.);
}
#[test]
fn spike_test() {
    let mut x1 = vec![10; 10];
    x1.push(1000);
    x1.append(&mut vec![10; 10]);
    let mut x2 = vec![10; 7];
    x2.push(1000);
    x2.append(&mut vec![10; 13]);
    let d = |x: &i32, y: &i32| (x - y).abs() as f32;
    let band = 15;
    let (score, _, _) = dtw_chiba(&x1, &x2, &d, band).unwrap();
    assert_eq!(score, 0.);
}

#[test]
fn short_band() {
    let mut x1: Vec<_> = (0..10).collect();
    x1.push(0);
    x1.push(0);
    let mut x2: Vec<_> = vec![0];
    x2.append(&mut (0..10).collect());
    x2.push(0);
    let d = |x: &i32, y: &i32| (x - y).abs() as f32;
    let band = 5;
    let (score, _, start) = dtw_chiba(&x1, &x2, &d, band).unwrap();
    assert_eq!(score, 0.);
    assert_eq!(start, 0);
}

#[test]
fn band_test() {
    let mut x1: Vec<_> = (0..1000).collect();
    x1.append(&mut vec![0; 10]);
    let mut x2 = vec![0; 9];
    x2.append(&mut (0..1000).collect());
    x2.push(0);
    let d = |x: &i32, y: &i32| (x - y).abs() as f32;
    let band = 21;
    let (score, _, start) = dtw_chiba(&x1, &x2, &d, band).unwrap();
    assert_eq!(score, 0.);
    assert_eq!(start, 0);
}
