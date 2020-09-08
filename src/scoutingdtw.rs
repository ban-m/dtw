use num::Num;
use num::Float;
use order_stat::kth_by;
const NUM_SCOUT:usize = 20;
const PACK_SIZE:usize = 5;
/// Scouting sub dynamic time warping with threshold bounding
pub fn scouting_threshold_dtw<D,F,T>(x1:&[D],x2:&[D],dist:&F,
                           num_scouts:Option<usize>,num_packs:Option<usize>,threshold:T) ->
    Result<(T,(Vec<usize>,Vec<usize>),usize),String>
    where F:Fn(&D,&D)->T, T:Float,D:Num+Copy,
{
    // x1 is query, x2 is reference.
    let x1len = x1.len();
    let num_packs = match num_packs{
        Some(res) => res,
        None => PACK_SIZE,
    };
    let num_scouts = match num_scouts{
        Some(res) => res,
        None => NUM_SCOUT,
    };
    let candidates = scouting(&x1[0..x1len/num_packs],x2,dist,num_scouts,threshold,x1len*3/2);
    if candidates.is_empty(){
        // early return.
        return Err("Early return".to_string())
    }
    let mut opt = T::infinity();
    let mut optcand = None;
    for (start,end) in candidates{
        if let Ok(res) = super::thresholddtw(x1,&x2[start..end],dist,threshold){
            if res.0 < threshold {
                return Ok(res)
            }else if res.0 < opt {
                opt = res.0;
                optcand = Some(res);
            }
        }
    }
    match optcand {
        Some(res) => Ok(res),
        None => Err("There's no candidates. This error should not happen.".to_string()),
    }
}

/// Scouting sub dynamic time warping.
/// This algorithm first compute "mini-" query to find candidates for "entire" query.
/// To determine start position, dynanic programming table is filled in a reverse order.
pub fn scouting_dtw<D,F,T>(x1:&[D],x2:&[D],dist:&F,
                           num_scouts:Option<usize>,num_packs:Option<usize>) ->
    Result<(T,(Vec<usize>,Vec<usize>),usize),String>
    where F:Fn(&D,&D)->T, T:Float,D:Num+Copy,
{
    scouting_threshold_dtw(x1,x2,dist,num_scouts,num_packs,T::infinity())
}

// Make candidates by mini query.
// After takeing top k position as "scouts",
// merging the region of two overlapping scouts,
// then sorting the scouts so that the main alignemnt can be
// executed from the most hopeful scout.
fn scouting<D,F,T>(x1:&[D],x2:&[D],dist:&F,num_scout:usize,threshold:T,subreflen:usize) -> Vec<(usize,usize)>
    where F:Fn(&D,&D)->T, T:Float,D:Num+Copy,
{
    let x1len = x1.len();// length of scout 
    let x2len = x2.len();
    if x1len == 0 || x2len == 0{
        //eprintln!("query:{},reference:{}",x1len,x2len);
        return vec![];
    }
    // make sure that this scouting procedure will be executed in
    // "reverse" manner.
    let bignum = T::epsilon().recip();
    let mut previous = vec![T::zero();x2len+1];
    let mut current = vec![bignum;x2len+1];
    let mut early_return;
    for i in 1..x1len + 1{
        early_return = true;
        for j in 1..x2len + 1{
            let d = dist(&x1[x1len-i],&x2[x2len-j]);// <- !!!!
            let min:T  = previous[j].min(previous[j-1].min(current[j-1]));
            current[j] = min + d;
            early_return = (current[j] > threshold) & early_return;
        }
        if early_return{
            //eprintln!("early return at :{}/{}",i,x1len);
            return vec![];
        }
        for k in 0..x2len+1{
            previous[k] = current[k];
            current[k] = bignum;// T::infinity();//
        }
    }
    // Make sure that the index is reversed order.
    // (score,index)
    let threshold = kth_by(&mut previous.clone(),num_scout,|x,y| x.partial_cmp(y).unwrap()).clone();
    let mut result:Vec<_> = previous.iter().enumerate()
        .filter_map(|(idx,&e)| if e <= threshold{Some((e,x2len-idx))}else{None})
        .collect();
    if result.is_empty(){
        vec![]
    }else{
        // sort the scouts in increasing order with respect to index;
        result.sort_by(|a,b|(a.1).cmp(&b.1));
        // remove overlapping region
        let mut result = remove_overlapping_scouts(&result,subreflen,x2len);
        // sort the scouts in increasing order with respect to score.
        result.sort_by(|a,b|(a.0).partial_cmp(&b.0).unwrap());
        result.into_iter().map(|a|(a.1,a.2)).collect()
    }
}

// (score,start,end)
fn remove_overlapping_scouts<T>(scouts:&Vec<(T,usize)>,subreflen:usize,maxlen:usize)-> Vec<(T,usize,usize)>
where T:Float{
    let mut res = vec![];
    let mut current_start_point = scouts[0].1;
    let mut current_end_point = (current_start_point + subreflen).min(maxlen);
    let mut current_best_score = scouts[0].0;
    for &scout in scouts{
        let start = scout.1;
        if start <= current_end_point{
            // scouts are overlapping. Merge them.
            current_best_score = (current_best_score).min(scout.0);
            current_end_point = (start + subreflen).min(maxlen);
        }else{
            // scuot is far enough from previous scout.
            // push the previous scout.
            assert!(current_start_point < current_end_point);
            res.push((current_best_score,current_start_point,current_end_point));
            current_start_point = if start < 20 {
                0
            }else{
                start - 20
            };
            current_end_point = (current_start_point + subreflen).min(maxlen);
            current_best_score = scout.0;
        }
    }
    res.push((current_best_score,current_start_point,current_end_point));
    res
}
