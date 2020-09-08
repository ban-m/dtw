use num::Float;
/// Wrapper struct invented by Keogh
#[derive(Debug)]
pub struct DynamicTimeWarping<D,T,F>
where F:Fn(&D,&D)->T, T:Float,D:Copy+PartialOrd{
    reference:Vec<D>,
    dist:F,
    lower_envelope:Vec<D>,
    upper_envelope:Vec<D>,
    bandwidth:usize,
    lowerbound_check_order:Vec<usize>,
}


impl<D,T,F> DynamicTimeWarping<D,T,F>
    where F:Fn(&D,&D)->T, T:Float,D:Copy+PartialOrd{
    /// constructor
    pub fn new(reference:Vec<D>,dist:F,bandwidth:usize)->DynamicTimeWarping<D,T,F>
    {
        let (lower_envelope,upper_envelope) = Self::envelope(&reference,bandwidth);
        let lowerbound_check_order = Self::ordering(&reference);
        DynamicTimeWarping{reference,
                           dist,
                           lower_envelope,
                           upper_envelope,
                           bandwidth,
                           lowerbound_check_order}
    }
    // compute envelope for given sequence
    fn envelope(events:&Vec<D>,bandwidth:usize)->(Vec<D>,Vec<D>){
        let len = events.len();
        let mut lower_envelope = vec![events[0];len];
        let mut upper_envelope = vec![events[0];len];
        for i in 0..events.len(){
            let start = if i < bandwidth { 0 } else{ i - bandwidth };
            let end = if i + bandwidth > len { len }else{i+bandwidth};
            for j in start .. end{
                if events[i] < lower_envelope[j] {
                    lower_envelope[j] = events[i];
                }
                if upper_envelope[j] < events[i]{
                    upper_envelope[j] = events[i];
                }
            }
        }
        (lower_envelope,upper_envelope)
    }
    // compute the order by which lower_bound_keogh() is executed.
    fn ordering(events:&Vec<D>)->Vec<usize>{
        let mut events:Vec<_> = events.iter().enumerate().collect();
        events.sort_by(|a,b|(b.1).partial_cmp(a.1).unwrap());
        events.iter().map(|e|e.0).collect()
    }
    /// Dynamic time warping for subsequence banded dynamic time warping.
    pub fn dtw(&self,query:&Vec<D>)->(T,usize){
        let mut best_so_far :T= T::infinity();
        let mut best_position = 0;
        let len = query.len();
        let (mut lb_q,mut lb_r,mut lb_m) = (0,0,0);
        let (query_lower_env,query_upper_env) = Self::envelope(query,self.bandwidth);
        let query_compare_ordering = Self::ordering(query);
        for (pos,subref) in self.reference.windows(query.len()).enumerate(){
            // no LB_Kim because its power is weak for very long query.
//            eprintln!("start position:{}",pos);
            let lb_query = self.lowerbound_keogh(subref,
                                                 &query_lower_env,
                                                 &query_upper_env,
                                                 &query_compare_ordering,
                                                 best_so_far);
//            eprintln!("lower bound:{:?}(query)",lb_query.to_f32());
            if lb_query > best_so_far{
                lb_q +=1;
                continue;
            }
            let lb_ref = self.lowerbound_keogh(query,
                                               &self.lower_envelope[pos..pos+len],
                                               &self.upper_envelope[pos..pos+len],
                                                &query_compare_ordering,
                                               best_so_far);
//            eprintln!("lower bound:{:?}(ref)",lb_ref.to_f32());
            if lb_ref > best_so_far{
                lb_r+=1;
                continue;
            }
            let cumulative_lb = if lb_ref < lb_query{
                self.cumulative_lower_bound(subref,
                                             &query_lower_env,
                                             &query_upper_env,
                                             lb_query)
            }else{
                self.cumulative_lower_bound(query,
                                             &self.lower_envelope[pos..pos+len],
                                             &self.upper_envelope[pos..pos+len],
                                             lb_ref)
            };
            let score = self.chiba_dtw_with_lower_bound(query,
                                                    subref,
                                                    &cumulative_lb,
                                                    best_so_far,
                                                    &mut lb_m);
            if score < best_so_far{
                best_so_far = score;
                best_position = pos;
            };
//            eprintln!("best so far {:?}",best_so_far.to_f32());
        }
        eprintln!("query:{},reference:{},while main dtw:{}",lb_q,lb_r,lb_m);
        (best_so_far,best_position)
    }
    fn lowerbound_keogh(&self,events:&[D],lower:&[D],upper:&[D],order:&[usize],
                        best_so_far:T)->T{
        let mut lb = T::zero();
        for &i in order{
//            eprint!("{},",i);
            if lb > best_so_far{
                return lb
            }else{
                lb = lb + if events[i] < lower[i]{
                    (self.dist)(&events[i],&lower[i])
                }else if events[i] > upper[i]{
                    (self.dist)(&events[i],&upper[i])
                }else{
                    T::zero()
                };
            }
        }
        lb
    }
    fn cumulative_lower_bound(&self,events:&[D],lower:&[D],upper:&[D],
                              lower_bound:T)->Vec<T>{
        // clb[i] is sum of the difference from i to the end.
        // First, compute lb, then compute clb[i] each.
        let mut lb = lower_bound;
        let mut clb = Vec::with_capacity(events.len());
        for i in 0..events.len(){
            lb = lb - if events[i] < lower[i]{
                (self.dist)(&events[i],&lower[i])
            }else if events[i] > upper[i]{
                (self.dist)(&events[i],&upper[i])
            }else{
                T::zero()
            };
            clb.push(lb);
        }
        clb
    }
    /// sub banded dynamic time warping using ucr optimization
    pub fn chiba_dtw_with_lower_bound(&self,query:&[D],subref:&[D],
                                  cumulative_lower_bound:&[T],best_so_far:T,
                                      lb_m:&mut u32)->T{
        // dynamic time warping. Not preserving traceback path.
        debug_assert!(query.len()==subref.len(),"r{},q{}",subref.len(),query.len());
//        eprintln!("{},{}",query.len(),subref.len());
        let len = query.len();
        let band = self.bandwidth;
        let mut previous = vec![T::infinity();2*band+1];
        let mut current = vec![T::infinity();2*band+1];
        //        eprintln!("{:?}\n{:?}",previous.len(),current.len());
        previous[band] = T::zero();
        for i in 1..len+1{
            let mut score = T::infinity();
            for j in 0..2*band+1{
                if i + j <= band || len + band <= i + j -1{
                    continue
                }else{
                    //eprint!("{},{}",i,j);
                    // gap transition on query
                    let gap = if j==0 {T::infinity()}else{current[j-1]};
                    // delete transition on query
                    let del = if j ==2*band {T::infinity()}else{previous[j+1]};
                    // match transision on query
                    let mat =  previous[j];
                    //eprintln!("-> score get");
                    current[j] = (self.dist)(&query[i-1],
                                             &subref[i+j-band-1]) +
                        if gap <= del && gap<= mat{
                            gap
                        }else if del <= gap && del <= mat{
                            del
                        }else{
                            mat
                        };
                    if current[j] < score{
                        score = current[j];
                    }
                }
            }
            if score+cumulative_lower_bound[i-1]>best_so_far{
                // early abandon
                *lb_m+=1;
                return score+cumulative_lower_bound[i-1]
            }
//            eprintln!();
            previous.clear();
            previous.append(&mut current);
            current.append(&mut vec![T::infinity();len]);
        }
        previous[band]
    }
}


