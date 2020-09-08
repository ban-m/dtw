use num::Float;

/// sub dynamic time warping with early abandoning.
/// This pruning idea comes from a simple observation that
/// once the score of dtw exceed a threshold, it
/// won't reach smaller score.
/// It returns Err(_) when it reaches threshold or
/// error.
#[inline]
pub fn thresholddtw<D,F,T>(x1:&[D],x2:&[D],dist:&F,threshold:T) 
                   -> Result<(T,(Vec<usize>,Vec<usize>),usize),String>
    where F:Fn(&D,&D)->T, T:Float
{
    // x1 is query,x2 is reference
    let x1len = x1.len();
    let x2len = x2.len();
    if x1len == 0 || x2len == 0{
        return Err("the length of given time series may be empty".to_string())
    }
    // dynamic programming path will be filled in
    // reference order.
    let bignum = T::epsilon().recip();
    let mut previous = vec![T::zero();x2len+1];
    let mut current = vec![bignum;x2len+1];
    let mut early_return;
    for i in 1..x1len + 1{
        early_return = true;
        for j in 1..x2len + 1{
            let min:T  = previous[j].min(previous[j-1].min(current[j-1]));
            current[j] = min+dist(&x1[i-1],&x2[j-1]);
            early_return = (current[j] > threshold) & early_return;
        }
        if early_return{
            return Err(format!("{},{}",i,x1len))
        }
        for k in 0..x2len+1{
            previous[k] = current[k];
            current[k] = bignum;
        }
    }
    Ok((previous.into_iter().fold(T::infinity(),|acc,x|if acc<x{acc}else{x}),
        (vec![],vec![]),0))
}

