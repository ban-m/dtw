use std::vec::Vec;
use std;
//const BIGNUM:f32 = 1000000.0;
use num::Float;
use num::Num;
use num;
use dtw_chiba;
use dtw_itakura;
use fastdtw;
use quickdtw;
use scoutingdtw;
/// Mode specifier to be used by other function 
/// to decide the dynamic time warping mode.
/// SakoeChiba(usize): Using Sakoe-Chiba band with the specified band width.
/// Itakura(usize): Using Itakura pentagram with the specified max band width.
#[derive(Debug,Copy,Clone)]
pub enum Mode where{
    /// Specifier for Sakoe-Chiba band.
    SakoeChiba(usize),
    /// Specifier for Itakuta's band.
    Itakura(usize),
    /// Specifier for sub dtw.
    Sub,
    /// Specifier for full dtw.
    Full,
    /// Specifier for fast dtw
    Fast(usize),
    /// Specifier for fast sub dtw
    FastSub(usize),
    /// Optimized sub dtw2
    QuickSub,
    /// Scouting sub dtw. Number of scouts and number of packs should be specified.
    Scouting(usize,usize),
}
impl std::fmt::Display for Mode{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result{
        let name = match *self{
            Mode::Fast(b) => format!("Fast{}",b),
            Mode::FastSub(b) => format!("FastSub,{}",b),
            Mode::Full => "Full".to_string(),
            Mode::Itakura(b) => format!("Itakura{}",b),
            Mode::QuickSub => "QuickSub".to_string(),
            Mode::SakoeChiba(b) => format!("SakoeChiba({})",b),
            Mode::Sub => "Sub".to_string(),
            Mode::Scouting(num_scouts,num_packs) => format!("Scouting({}_{})",num_scouts,num_packs),
        };
        write!(f, "{}", name)
    }
}


/// A generic function of dymanic time warping.
/// This function can execute dynamic time warping 
/// for almost all data type with appropriate function dist.
/// (score,(path,path),index to start the alignment)
/// X1 IS QUERY AND X2 IS REFERENCE NOT VISE VERSA
pub fn dtw<D,F,T>(x1:&[D],x2:&[D],mode:Mode,dist:&F) -> Result<(T,(Vec<usize>,Vec<usize>),usize),String>
    where F:Fn(&D,&D)->T, T:Float,D:Num+Copy,
{
    let err = format!("invalid input: the lengths are different,while restricted mode was chosen:{},{}",x1.len(),x2.len()).to_string();
    match mode{
        Mode::SakoeChiba(_) if x1.len() != x2.len() => Err(err),
        Mode::SakoeChiba(n) => dtw_chiba::dtw_chiba(x1,x2,dist,n),
        Mode::Itakura(_) if x1.len() != x2.len() => Err(err),
        Mode::Itakura(n) => dtw_itakura::dtw_itakura(x1,x2,dist,n),
        Mode::Full => dtw_norestrict(x1,x2,dist,mode),
        Mode::Sub  => dtw_norestrict(x1,x2,dist,mode),
        Mode::Fast(radius) => fastdtw::fast_dtw(x1,x2,dist,radius,false),
        Mode::FastSub(radius) => fastdtw::fast_dtw(x1,x2,dist,radius,true),
        Mode::QuickSub => quickdtw::quickdtw(x1,x2,dist),
        Mode::Scouting(num_scouts,num_packs) => scoutingdtw::scouting_dtw(x1,x2,dist,Some(num_scouts),Some(num_packs)),
    }
}

fn dtw_norestrict<D,F,T>(x1:&[D],x2:&[D],dist:&F,mode:Mode) 
                   -> Result<(T,(Vec<usize>,Vec<usize>),usize),String> 
    where F:Fn(&D,&D)->T, T:Float
{
    // x1 is qyery,x2 is reference
    let x1len = x1.len();
    let x2len = x2.len();
    if x1len == 0 || x2len == 0{
        return Err("the length of given time series may be empty".to_string())
    }
    let mut res = match mode {
        Mode::Full => vec![vec![T::infinity();x2len+1];x1len+1],
        Mode::Sub => { let mut dp = vec![vec![T::infinity();x2len+1];x1len+1];
                       for j in 0..x2len {
                           dp[0][j] = num::zero();
                       };
                       dp},
        _ => unreachable!(),
    };
    res[0][0] = num::zero();
    for i in 1..x1len + 1{
        for j in 1..x2len + 1{
            let cost = dist(&x1[i-1],&x2[j-1]);
            let min:T  = res[i][j-1].min(res[i-1][j].min(res[i-1][j-1]));
            res[i][j] = if !min.is_infinite(){ min + cost } else { min };
        }
    }
    let mut x1path = vec![];
    let mut x2path = vec![];
    let mut x1pos = x1len;
    let (mut x2pos,score) = match mode {
        Mode::Full => (x2len,res[x1len][x2len]),
        Mode::Sub => res[x1len].iter()
            .enumerate()
            .fold((0,Float::infinity()),|(idx,acc),(i,&score)| if acc < score {
                (idx,acc)
            }else{
                (i,score)}),
        _ => unreachable!(),
    };
    while x1pos > 0 {
        x1path.push(x1pos-1);
        x2path.push(x2pos-1);
        let next = nextpos(x1,x2,&res,x1pos,x2pos,dist);
        x1pos = next.0;
        x2pos = next.1;
    }
    x1path.reverse();
    x2path.reverse();
    match mode {
        Mode::Full => Ok((score,(x1path,x2path),0)),
        Mode::Sub => Ok((score,(x1path,x2path),x2pos-1)),
        _ => unreachable!("could not happen"),
    }
}

fn nextpos<D,F,T>(x1:&[D],
                  x2:&[D],
                  res:&Vec<Vec<T>>,
                  x1pos:usize,
                  x2pos:usize,dist:&F) ->  (usize,usize) 
    where F:Fn(&D,&D)->T, T:Float{
    let current_score = res[x1pos][x2pos] - dist(&x1[x1pos-1],&x2[x2pos-1]);
    let match_path = res[x1pos-1][x2pos-1] - current_score;
    let gap_path = res[x1pos][x2pos-1] - current_score;
    let del_path = res[x1pos-1][x2pos] - current_score;
    if match_path < gap_path && match_path < del_path{
        (x1pos-1,x2pos-1)
    }else if gap_path < match_path && gap_path < del_path {
        (x1pos,x2pos-1)
    }else{
        (x1pos-1,x2pos)
    }
}


// fn dist(x:f32,y:f32)->f32{
//     (x-y).powi(2)
// }

/// Dynamic time warping for 1 dimentional data.
/// 
/// Because the naive PD algorithm is used,
/// the time complexity is O(mn) where m and n are the lengths 
/// of given data.
/// Currently no "window" mode is implemented.
pub fn dtw_with_path(x1:&[f32],x2:&[f32]) -> (f32,(Vec<usize>,Vec<usize>)){
    let d = |x:&f32,y:&f32| (x-y).powi(2);
    let (score,(x1path,x2path),_) = dtw(x1,x2,Mode::Full,&d).unwrap();
    // res[0][0] = 0.0;
    // for i in 1..x1len+1{ 
    //     for j in 1..x2len+1{
    //         let cost = dist(x1[i-1],x2[j-1]);
    //         res[i][j] = cost + res[i-1][j].min(res[i][j-1].min(res[i-1][j-1]));
    //     }
    // }
    //    (res[x1len][x2len],res)
    (score,(x1path,x2path))
}

/// Dynamic time warping for 1 dimentional data.
///
/// This is a wrapper function of dtw_with_path currently.
/// Maybe this function will be modified so that 
/// the time complexity is strictly less than O(mn).
pub fn normal_dtw(x1:&[f32],x2:&[f32]) -> f32{
    dtw_with_path(x1,x2).0
}

// ///Dynamic time warping for 2 dimentional data.
// ///
// ///A naive implementation for dynamic time warping, thus 
// ///the time complexity is O(mn),where m and n are the lengths of the data respectively.
// ///If you only want to compute the dtw score, use dtw(x1,x2).
// pub fn dtw2_with_path(x1:&[[f32;2]],x2:&[[f32;2]]) -> (f32,(Vec<usize>,Vec<usize>)){
//     let d = |x:&[f32;2],y:&[f32;2]| ((x[0]-y[0]).powi(2) + (x[1]-y[1]).powi(2)).sqrt();
//     let (score,(x1path,x2path),_) = dtw(x1,x2,Mode::Full,&d).unwrap();
//     (score,(x1path,x2path))
// }

// /// Dynamic time warping for 2 dimentional data.
// ///
// /// This is a wrapper function of dtw2_with_path.
// /// Almost the same as dtw().
// ///
// pub fn dtw2(x1:&[[f32;2]],x2:&[[f32;2]]) -> f32{
//     dtw2_with_path(x1,x2).0
// }

/// Subsequence dynamic time warping for 1 dimentional data.
/// x1 is query,x2 is reference
pub fn subdtw_with_path(x1:&[f32],x2:&[f32])->(f32,(Vec<usize>,Vec<usize>)){
    let d = |x:&f32,y:&f32| (x-y).powi(2);
    let (score,(x1path,x2path),_) = dtw(x1,x2,Mode::Sub,&d).unwrap();
    (score,(x1path,x2path))
    // let x1len = x1.len();
    // let x2len = x2.len();
    // let mut res = vec![vec![BIGNUM;x2len+1];x1len+1];
    // for j in 0..x2len+1{
    //     res[0][j] = 0.0;
    // }
    // for i in 1..x1len+1{
    //     for j in 1..x2len+1{
    //         let cost = dist(x1[i-1],x2[j-1]);
    //         res[i][j] = cost + res[i-1][j].min(res[i][j-1].min(res[i-1][j-1]));
    //     }
    // }
    // (res[x1len].iter().fold(BIGNUM,|acc,&x| acc.min(x)),res)
}


/// Subsequnece dynamic time warping for 1 dimentional data.
///
/// First of all, the first argument is query whereas x2 is reference. In other words, x1 is global,x2 is local.
/// Note that the naive DP algorithm is used and the time complexity 
/// is not changed at all: still O(mn).
/// If you want to know where to start the alignment, use subdtw_map() instead.
pub fn subdtw(x1:&[f32],x2:&[f32])->f32{
    subdtw_with_path(x1,x2).0
}

// fn nextposition(x1:&[f32],x2:&[f32],score:&Vec<Vec<f32>>,x1pos:usize,x2pos:usize)->(usize,usize){
//     let current_score = score[x1pos][x2pos] - dist(x1[x1pos-1],x2[x2pos-1]);
//     let match_path = score[x1pos-1][x2pos-1] - current_score;
//     let gap_path = score[x1pos][x2pos-1] - current_score;
//     let del_path = score[x1pos-1][x2pos] - current_score;
//     if match_path < gap_path && match_path < del_path {
//         (x1pos-1,x2pos-1)
//     }else if gap_path < match_path && gap_path < del_path{
//         (x1pos,x2pos-1)
//     }else{
//         (x1pos-1,x2pos)
//     }
// }

/// Compute the subsequence dynamic time warping path from score matrix
///
/// Note that the collumn should correspond to query, the row is the other(reference).
/// To get the warp itself, you need to access the argument after calling 
/// this function like: x1[res[0]],x1[res[1]],...
///
pub fn sub_warppath(x1:&[f32],x2:&[f32],_score:&Vec<Vec<f32>>)-> (Vec<usize>,Vec<usize>){
    // find the starting position
    let d = |x:&f32,y:&f32| (x-y).powi(2);
    let (_,(x1path,x2path),_) = dtw(x1,x2,Mode::Full,&d).unwrap();
    (x1path,x2path)
    //     let mut query_pos = x1.len();
    // let mut ref_pos = score[x1.len()].iter()
    //     .enumerate()
    //     .fold((BIGNUM,0),|(acc,idx),(i,&x)| if acc < x { (acc,idx)}else{(x,i)}).1;
    // let end = 0;
    // let mut querypath = vec![];
    // let mut referencepath = vec![];
    // while query_pos > end {
    //     assert!(ref_pos > 1);
    //     querypath.push(query_pos-1);
    //     referencepath.push(ref_pos-1);
    //     let next_pos:(usize,usize) = nextposition(x1,x2,score,query_pos,ref_pos);
    //     query_pos = next_pos.0;
    //     ref_pos = next_pos.1;
    // }
    // querypath.reverse();
    // referencepath.reverse();
    // (querypath,referencepath)
}
#[allow(dead_code)]
fn dp_flush<T>(dp:&Vec<Vec<T>>) where
    T:Float
{
    let row = dp.len();
    let column = dp[0].len();
    for i in 0..row{
        for j in 0..column{
            if dp[i][j].is_infinite(){
                print!("inf");
            }else{
                print!("{:3.0}",dp[i][j].to_f32().unwrap());
            }
        }
        println!("");
    }
}

/// Subsequence dynamic time warping for 1 dimentional data.
///
/// This function returns not only the optimal score,
/// but also the reference index from which the alignment start.
/// Although the simple DP algorithm is used, there is 
/// no 'trackbacking', so the compute time is almost the same
/// as other method such as subdtw().
pub fn subdtw_map(x1:&[f32],x2:&[f32])->(f32,usize){
    let d = |x:&f32,y:&f32| (x-y).powi(2);
    let (score,_,location) = dtw(x1,x2,Mode::Sub,&d).unwrap();
    (score,location)
//     let x1len = x1.len();
//     let x2len = x2.len();
// //    let mut res = vec![vec![BIGNUM;x2len+1];x1len+1];
//     let mut res = vec![BIGNUM;(x2len+1)*(x1len+1)];
//     for j in 0..x2len+1{
//         res[j] = 0.0;
//     }
//     for i in 1..x1len+1{
//         for j in 1..x2len+1{
//             let cost = dist(x1[x1len-i],x2[x2len-j]);
//             res[i*(x2len+1)+j] = cost + res[(i-1)*(x2len+1)+j].min(res[i*(x2len+1)+j-1].min(res[(i-1)*(x2len+1)+j-1]));
//         }
//     }
//     res.iter().skip(x1len*(x2len+1)).enumerate().fold((BIGNUM,0),|(min,minidx),(idx,&x)|{
//         if x < min {
//             (x,idx)
//         }else{
//             (min,minidx)
//         }
//     })
}

// /// Subsequence dynamic time warping for 1 dim data by specified distance measure
// pub fn subdtw_by<F>(x1:&[f32],x2:&[f32],dist:F) -> (f32,Vec<Vec<f32>>)
//     where F : Fn(f32,f32) -> f32 {
//      let x1len = x1.len();
//     let x2len = x2.len();
//     let mut res = vec![vec![BIGNUM;x2len+1];x1len+1];
//     for j in 0..x2len+1{
//         res[0][j] = 0.0;
//     }
//     for i in 1..x1len+1{
//         for j in 1..x2len+1{
//             let cost = dist(x1[i-1],x2[j-1]);
//             res[i][j] = cost + res[i-1][j].min(res[i][j-1].min(res[i-1][j-1]));
//         }
//     }
//     (res[x1len].iter().fold(BIGNUM,|acc,&x| acc.min(x)),res)
// }
