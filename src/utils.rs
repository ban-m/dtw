use std;
#[derive(Debug,Clone,Copy)]
/// This enum is for preprocess of dynamic time wapring.
/// When optimal_dtw is executed with Flat enum,
/// histogram equilization for query is done before 
/// main algorithm carries out.
pub enum Prep{
    /// For histogram equilization. Make sure only query will be equilized and
    /// to equilization cumulative distribution function is needed.
    Flat,
    /// Do nothing for preprocess. Usually one want to use this.
    Normal,
}

impl Prep{
    /// "FLAT" or "NORMAL" can be used to make a instance.
    pub fn new(prep:&str) -> Option<Prep>{
        match prep{
            "FLAT" => Some(Prep::Flat),
            "NORMAL" => Some(Prep::Normal),
            _ => None,
        }
    }
}

impl std::fmt::Display for Prep{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result{
        let name = match self{
            &Prep::Flat => "Flat",
            &Prep::Normal => "Normal",
        };
        write!(f, "{}", name)
    }
}

#[derive(Debug,Clone,Copy)]
/// Enum for execute optimal dynamic time wapring method.
/// Currently, four methods can be determined.
pub enum Method{
    /// Chiba Sakoe band with normal squared-Euclidian metric(bandwidth = query size / 10).
    ChibaNormal,
    /// Chiba Sakoe band with hill function metric (bandwidth = querysize/10).
    ChibaHill,
    /// Sub dynamic time warping with normal squared-Euclidan metric(bandwidth = querysize/10)
    SubNormal,
    /// Sub dybamic time warping with hill function metric(bandwidth = querysize/10)
    SubHill,
}

///parse mode. Currently parse only Sub, Chiba and Scouting dtw.
pub fn get_mode(mode:&str)->std::result::Result<super::Mode,String>{
    if mode.starts_with("Sub"){
        Ok(super::Mode::Sub)
    }else if mode.starts_with("Chiba"){
        let bandwidth:usize = match mode.split(',').nth(1)
            .and_then(|e|e.parse().ok()){
                Some(b) => b,
                None => return Err("Given chiba but couldn't parse bandwidth correctly"
                                   .to_string()),
            };
        Ok(super::Mode::SakoeChiba(bandwidth))
    }else if mode.starts_with("Scouting"){
        let contents:Vec<usize> = mode.split(',').skip(1)
            .filter_map(|e|e.parse().ok()).collect();
        if contents.len() == 2 {
            Ok(super::Mode::Scouting(contents[0],contents[1]))
        }else{
            return Err("Given scouting but couldn't parse argument.".to_string())
        }
    }else if mode.starts_with("Fast"){ 
        let bandwidth:usize = match  mode.split(",")
            .nth(1)
            .and_then(|e| e.parse().ok()){
                Some(b) => b,
                None => return Err(format!("Given chiba but couldn't parse bandwidth correctly{}",mode))
            };
        Ok(super::Mode::Fast(bandwidth))
    }else{
        return Err(format!("invalid mode name:{}",mode))
    }
}


impl Method{
    /// Constructor of Method enum.
    pub fn new(method:&str)->Option<Method>{
        match method{
            "CHIBANORMAL" => Some(Method::ChibaNormal),
            "CHIBAHILL" => Some(Method::ChibaHill),
            "SUBNORMAL" => Some(Method::SubNormal),
            "SUBHILL" => Some(Method::SubHill),
            _ => None,
        }
    }
}
impl std::fmt::Display for Method{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result{
        let name = match self{
            &Method::ChibaHill => "ChibaHill",
            &Method::ChibaNormal => "ChibaNormal",
            &Method::SubHill => "SubHill",
            &Method::SubNormal => "SubNormal",
        };
        write!(f, "{}", name)
    }
}


#[inline]
fn hill(x:&f32,y:&f32) -> f32{
    let d = (x-y).powi(2);
    d/(1. + d)
}

#[inline]
fn normal(x:&f32,y:&f32) -> f32{
    (x-y).powi(2)
}

#[inline]
fn padding(e:&f32)->Vec<f32>{
    use rand::thread_rng;
    use rand::Rng;
    let mut res = vec![e.clone()];
    let continue_prob = 0.84/1.84;
    let mut rng = thread_rng();
    while rng.gen_range(0.,1.) < continue_prob{
        res.push(e.clone());
    }
    res
}


#[inline]
fn padding_reference(reference:&[f32]) -> Vec<f32>{
    reference.iter().map(|e| padding(e)).fold(vec![],|mut acc,mut x|{acc.append(&mut x);acc})
}

/// fast dtw with padding
pub fn padding_fast_dtw<F>(query:&[f32],reference:&[f32],dist:&F,raidus:usize,is_sub:bool)
                           ->Result<(f32,(Vec<usize>,Vec<usize>),usize),String>
    where F:Fn(&f32,&f32)->f32{
    let reference = padding_reference(reference);
    super::fast_dtw(query,&reference,dist,raidus,is_sub)
}

/// optimal dynamic time warping.
/// Make sure if you specified prep as Prep::Flat, you should 
/// super::histogram_modify() to reference sequence before and
/// give this function its cumulative distribution function.
pub fn optimal_dtw(query:&Vec<f32>,reference:&Vec<f32>,method:&Method,
                   prep:&Prep,cdf:&Vec<f32>)->Option<f32>{
    const BANDWIDTH:usize = 51;
    let query =  match prep {
        &Prep::Flat => super::histgram_modify(&query,cdf),
        &Prep::Normal => query.clone(),
    };
    match method{
        &Method::SubHill => super::dtw(&query,&reference,
                                       super::Mode::Sub,&hill)
            .map(|e|e.0).ok(),
        &Method::SubNormal => super::dtw(&query,&reference,
                                         super::Mode::Sub,&normal)
            .map(|e|e.0).ok(),
        &Method::ChibaHill | &Method::ChibaNormal =>{
            let mut skip_count = 0;
            let mut best_so_far = 10000.;
            let is_hill = match method {
                &Method::ChibaHill => true,
                _ => false};
            for subref in reference.windows(query.len()){
                if skip_count > 0{
                    skip_count -=1;
                    continue;
                }
                let res = if is_hill{
                    super::dtw(&query,subref,super::Mode::SakoeChiba(BANDWIDTH),&hill)
                }else{
                    super::dtw(&query,subref,super::Mode::SakoeChiba(BANDWIDTH),&normal)
                };
                if let Ok(res) = res{
                    if res.0 < best_so_far{
                        best_so_far = res.0;
                        skip_count += BANDWIDTH/2;
                    }else{
                        skip_count += query.len()/10;
                    }
                }
            }
            Some(best_so_far)
        }
    }
}


/// optimal dynamic time warping.
/// Make sure if you specified prep as Prep::Flat, you should 
/// super::histogram_modify() to reference sequence before and
/// give this function its cumulative distribution function.
pub fn folded_optimal_dtw(query:&[f32],reference:&[f32],method:&Method,
                   prep:&Prep,cdf:&Vec<f32>)->Option<f32>{
    let subquery:Vec<_> = query.chunks(2)
        .map(|e| if e.len() == 2{(e[0]+e[1])/2.}else{e[0]})
        .collect();
    let subref = reference.chunks(2)
        .map(|e| if e.len() == 2{(e[0]+e[1])/2.}else{e[0]})
        .collect();
    let (_,location) = match optimal_dtw_with_location(&subquery,&subref,method,prep,cdf){
        Some(res) => res,
        None => return None,
    };
    let query =  match prep {
        &Prep::Flat => super::histgram_modify(&query.to_vec(),cdf),
        &Prep::Normal => query.to_vec(),
    };
    let querysize = query.len();
    let refsize = reference.len();
    let start = if location < querysize{0}else{location-querysize};
    let end = if location + querysize > refsize{
        refsize
    }else{
        location+querysize
    };
    let score = match method{
        &Method::SubHill => super::dtw(&query,&reference[start..end],
                                       super::Mode::Sub,&hill).unwrap().0,
        &Method::SubNormal => super::dtw(&query,&reference[start..end],
                                         super::Mode::Sub,&normal).unwrap().0,
        &Method::ChibaHill => chiba_optimal_dtw(&query,
                                                &reference[start..end],
                                                &super::Mode::SakoeChiba(51),true)
            .unwrap(),
        &Method::ChibaNormal =>chiba_optimal_dtw(&query,
                                                 &reference[start..end],
                                                 &super::Mode::SakoeChiba(51),false)
            .unwrap(),
    };
    Some(score)
}

fn optimal_dtw_with_location(query:&Vec<f32>,reference:&Vec<f32>,method:&Method,
                                 prep:&Prep,cdf:&Vec<f32>)->Option<(f32,usize)>{
    let query =  match prep {
        &Prep::Flat => super::histgram_modify(&query,cdf),
        &Prep::Normal => query.clone(),
    };
    let (score,_,location) = match method{
        &Method::SubHill => super::dtw(&query,&reference,
                                       super::Mode::Sub,&hill).unwrap(),
        &Method::SubNormal => super::dtw(&query,&reference,
                                         super::Mode::Sub,&normal).unwrap(),
        &Method::ChibaHill | &Method::ChibaNormal =>{
            let querysize = query.len();
            let refsize = reference.len();
            let mut result = None;
            let mut opt = 1000000.;
            let mut offset = 0;
            let bandwidth = 51;
            while offset+querysize < refsize {
                let subref = &reference[offset..offset+querysize];
                let subref = &padding_reference(subref)[0..querysize];
                let (score,t,sub_start) = match method {
                    &Method::ChibaHill => super::dtw(&query,
                                                     subref,
                                                     super::Mode::SakoeChiba(bandwidth),
                                                     &hill).unwrap(),
                    &Method::ChibaNormal => super::dtw(&query,
                                                       subref,
                                                       super::Mode::SakoeChiba(bandwidth),
                                                       &normal).unwrap(),
                    _ => unreachable!(),
                };
                if opt > score {
                    result = Some((score,t,offset+sub_start));
                    opt = score;
                    offset += bandwidth/2;
                }else{
                    offset += querysize/10;
                }
            };
            match result{
                Some(res) => res,
                None => unreachable!("no!!!"),
            }
        }
    };
    Some((score,location))
}




///fast mode
pub fn chiba_skipping_dtw(query:&[f32],reference:&[f32],
                                  mode:super::Mode,is_hill:bool)->Option<f32>{
    let mut best_so_far = 100000.;
    let mut skip_count = 0;
    let querysize = query.len();
    let bandwidth = match mode{
        super::Mode::SakoeChiba(b) => b,
        _ => {eprintln!("not valid input");return None},
    };
    
    for subref in reference.windows(query.len()){
        if skip_count > 0{
            skip_count -=1;
            continue;
        }
        let res = if is_hill{
            super::dtw(query,subref,mode,&hill)
        }else{
            super::dtw(query,subref,mode,&normal)
        };
        if let Ok(res) = res{
            if res.0 < best_so_far{
                best_so_far = res.0;
                skip_count = bandwidth/2;
            }else{
                skip_count = querysize/10;//2*bandwidth;
            }
        }
    }
    Some(best_so_far)
}

/// Debug mode,too slow
pub fn chiba_optimal_dtw(query:&[f32],reference:&[f32],
                         mode:&super::Mode,is_hill:bool)->Option<f32>{
    let querysize = query.len();
    let bandwidth = match mode{
        &super::Mode::SakoeChiba(b) => b,
        _ => {eprintln!("not valid input");return None},
    };
    let mut opt = 1000000.;
    for subref in reference.windows(querysize){
        let subref = &padding_reference(&subref)[0..querysize];
        let result = if is_hill{
            super::dtw(&query,
                       subref,
                       super::Mode::SakoeChiba(bandwidth),
                       &hill)
        }else{
            super::dtw(&query,
                       subref,
                       super::Mode::SakoeChiba(bandwidth),
                       &normal)
        };
        if let Ok((score,_,_)) = result{
            if opt > score {
                opt = score;
            }
        }
    }
    Some(opt)
}

///optimal dtw. Note this function does not need any "util" enum.
pub fn dtw_wrapper(query:&[f32],reference:&[f32],mode:&super::Mode,metric:&str,
                   prep:&Option<Vec<f32>>,threshold:&Option<f32>)->Option<f32>{
    use super::Mode;
    let query = match prep {
        &Some(ref cdf) => super::histgram_modify(query,cdf),
        &None => query.to_vec(),
    };
    let is_hill = metric == "Hill" || metric == "hill" || metric == "HILL";
    match mode{
        &Mode::Sub | &Mode::QuickSub => match threshold{
            &Some(threshold) => if is_hill{
                super::thresholddtw(&query,&reference,
                                    &hill,threshold)
            }else{
                super::thresholddtw(&query,&reference,
                                    &normal,threshold)
            },
            &None => if is_hill{
                super::dtw(&query,&reference,mode.clone(),&hill)
            }else{
                super::dtw(&query,&reference,mode.clone(),&normal)
            },
        }.map(|e|e.0).ok(),
        &Mode::Scouting(scouts,packs) => match threshold{
            &Some(threshold) => if is_hill{
                super::scouting_threshold_dtw(&query,&reference,
                                              &hill,
                                              Some(scouts),Some(packs),threshold)
            }else{
                super::scouting_threshold_dtw(&query,&reference,
                                              &normal,
                                              Some(scouts),Some(packs),threshold)
            },
            &None => if is_hill{
                super::dtw(&query,&reference,mode.clone(),&hill)
            }else{
                super::dtw(&query,&reference,mode.clone(),&normal)
            }
        }.map(|e|e.0).ok(),
        &Mode::SakoeChiba(b) =>{
            let mut skip_count = 0;
            let mut best_so_far = 1000000.;
            for subref in reference.windows(query.len()){
                if skip_count > 0{
                    skip_count -=1;
                    continue;
                }
                let res = if is_hill{
                    super::dtw(&query,subref,mode.clone(),&hill)
                }else{
                    super::dtw(&query,subref,mode.clone(),&normal)
                };
                if let Ok(res) = res{
                    if res.0 < best_so_far{
                        best_so_far = res.0;
                        skip_count += b;
                    }else{
                        skip_count += query.len()/10;
                    }
                }
            }
            Some(best_so_far)
        }
        &Mode::FastSub(r) => super::fast_dtw(&query,&reference,&normal,r,true)
            .map(|e|e.0).ok(),
        _ => None,
    }
    
}


