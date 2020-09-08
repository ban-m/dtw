use num::Float;

/// quick DTW. Sub dynamic time warping
pub fn quickdtw<D,F,T>(x1:&[D],x2:&[D],dist:&F) 
                   -> Result<(T,(Vec<usize>,Vec<usize>),usize),String> 
    where F:Fn(&D,&D)->T, T:Float
{

    let x1len = x1.len();
    let x2len = x2.len();
    if x1len == 0 || x2len == 0{
        return Err("the length of given time series may be empty".to_string())
    }
    let zero = T::zero();
    let inf = T::infinity();
    let mut previous = vec![inf;x1len+1];
    let mut current = vec![inf;x1len+1];
    previous[0] = zero;
    current[0] = zero;
    let mut opt_score = inf;
    for j in 1..x2len + 1{
        for i in 1..x1len + 1{
            let cost = dist(&x1[i-1],&x2[j-1]);
            let min:T  = previous[i].min(previous[i-1].min(current[i-1]));
            current[i] = if !min.is_infinite(){ min + cost } else { min };
        }
        opt_score = opt_score.min(current[x1len]);
        previous.clear();
        previous.append(&mut current);
        current = vec![inf;x1len+1];
        current[0] = zero;
    }
    Ok((opt_score,(vec![],vec![]),0))
}
