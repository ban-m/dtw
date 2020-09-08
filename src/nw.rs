use num::Float;

/// NW algorithm interpreted by dtw.
pub fn needleman_dtw<D,F,T>(x1:&[D],x2:&[D],dist:&F,gap:T) 
                   -> Result<(T,(Vec<usize>,Vec<usize>),usize),String> 
    where F:Fn(&D,&D)->T, T:Float
{

    let x1len = x1.len();
    let x2len = x2.len();
    if x1len == 0 || x2len == 0{
        return Err("the length of given time series may be empty".to_string())
    }
    let mut previous = vec![T::infinity();x1len+1];
    let mut current = vec![T::infinity();x1len+1];
    previous[0] = T::zero();
    current[0] = T::zero();
    let mut opt_score = T::infinity();
    for j in 1..x2len + 1{
        for i in 1..x1len + 1{
            let cost = dist(&x1[i-1],&x2[j-1]);
            current[i]  = (previous[i]+gap).min((previous[i-1]+cost).min(current[i-1]+gap));
        }
        opt_score = opt_score.min(current[x1len]);
        previous.clear();
        previous.append(&mut current);
        current = vec![T::infinity();x1len+1];
        current[0] = T::zero();
    }
    Ok((opt_score,(vec![],vec![]),0))
}
