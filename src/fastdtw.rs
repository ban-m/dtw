
use num::Float;
use num::Num;
use std::collections::HashMap;
/// O(N) implementation of dynamic time warping.
///
///
pub fn fast_dtw<D,F,T>(x1:&[D],x2:&[D],dist:&F,radius:usize,is_sub:bool)->
    Result<(T,(Vec<usize>,Vec<usize>),usize),String>
    where F:Fn(&D,&D) -> T,T:Float,D:Num+Copy
// when is_sub is true, x1 is query and x2 is reference.
{
    let min_length = radius + 2;
    if x1.len() < min_length || x2.len() < min_length{
        // normal dtw.
        // full window
        let window = (0..x1.len()).map(|i| (0..x2.len()).map(|j|(i+1,j+1)).collect())
            .fold(Vec::with_capacity(min_length*min_length),|mut acc,mut x|{acc.append(&mut x);acc});
        let res = window_dtw(x1,x2,dist,&window,is_sub);
        res
    }else{
        // recursive call
        let x1_shrinked = reduce_by_half(x1);
        let x2_shrinked = reduce_by_half(x2);
        let (_score,(x1path,x2path),_location) = fast_dtw(&x1_shrinked,&x2_shrinked,dist,radius,is_sub)?;
        let window = expand_window(&x1path,&x2path,x2.len(),x1.len(),radius);
        window_dtw(x1,x2,dist,&window,is_sub)
    }
}

fn window_dtw<D,F,T>(x1:&[D],x2:&[D],dist:&F,window:&[(usize,usize)],is_sub:bool) -> 
    Result<(T,(Vec<usize>,Vec<usize>),usize),String>
    where F:Fn(&D,&D) -> T,T:Float,D:Num+Copy{
    let mut dp = HashMap::new();//dp table. This map ontains 1:optimal score,2:traceback information
    for &(i,j) in window{
//        eprint!("({},{})->",i,j);
        let d = dist(&x1[i-1],&x2[j-1]);
        let paths = get_previous_scores(&dp,i,j,is_sub);
        let opt = get_min(paths,d);
        dp.insert((i,j),opt);
//        eprintln!("{:?}",opt.0.to_f32());
    }
    let (opt,max_position) = if is_sub{
        get_optimal(&dp,&window,x1.len())
    }else{
        (dp.get(&(x1.len(),x2.len())).ok_or("error occured while extracting optimal score".to_string())?.0,x2.len())
    };
    let (start_position,x1path,x2path) = get_paths(&dp,x1.len(),max_position);
    // for (&i,&j) in x1path.iter().zip(x2path.iter()){
    //     eprint!("({},{})->",i,j);
    // }
    // eprintln!("finish");
    Ok((opt,(x1path,x2path),start_position))
}

// #[inline]
// fn print_path(p1:&[usize],p2:&[usize]){
//     p1.iter().zip(p2.iter()).map(|(&i,&j)| eprint!("({},{})->",i,j)).count();
//     eprintln!();
// }

#[inline]
fn get_previous_scores<F>(dp:&HashMap<(usize,usize),(F,usize,usize)>,i:usize,j:usize,is_sub:bool)->
    [(F,usize,usize);3] where F:Float
{
    let p1 = get_score(&dp,i-1,j-1,is_sub);// "match" path
    let p2 = get_score(&dp,i-1,j,is_sub);//"insert" path
    let p3 = get_score(&dp,i,j-1,is_sub);//"gap" path
    [(p1,i-1,j-1),(p2,i-1,j),(p3,i,j-1)]
}

#[inline]
fn get_score<F>(dp:&HashMap<(usize,usize),(F,usize,usize)>,i:usize,j:usize,is_sub:bool)->F where F:Float{
    // get the optimal scores of previous cells
    match dp.get(&(i,j)){
        None if i == 0 && j == 0 && !is_sub => F::zero(),// initial score of full dtw
        None if i == 0 && is_sub => F::zero(), // initial score of sub dtw
        None => Float::infinity(),
        Some(entry) => entry.0
    }
}

#[inline]
fn get_min<F>(paths:[(F,usize,usize);3],prev:F)->(F,usize,usize) where F:Float{
    let opt = if paths[0].0 <= paths[1].0 && paths[0].0 <= paths[2].0{
        paths[0]
    }else if paths[1].0 <= paths[0].0 && paths[1].0 <= paths[2].0{
        paths[1]
    }else{
        paths[2]
    };
    if opt.0.is_infinite(){
        opt
    }else{
        (opt.0+prev,opt.1,opt.2)
    }
}

#[inline]
fn get_optimal<F>(dp:&HashMap<(usize,usize),(F,usize,usize)>,window:&[(usize,usize)],x1len:usize)->(F,usize)
    where F:Float
{
    window.iter().filter(|&&(i,_)| i==x1len)
        .filter_map(|&(i,j)| dp.get(&(i,j)).map(|&(opt,_,_)|(opt,j)))
        .fold((Float::infinity(),0),|acc,x| if acc.0 < x.0 { acc } else { x })
}

#[inline]
fn get_paths<F>(dp:&HashMap<(usize,usize),(F,usize,usize)>,i:usize,j:usize)->(usize,Vec<usize>,Vec<usize>)
    where F:Float
{
    let (mut x1path,mut x2path)  = (vec![],vec![]);
    let (mut i,mut j) = (i,j);
    while let Some(&(_,x1,x2)) = dp.get(&(i,j)){
        x1path.push(x1);
        x2path.push(x2);
        i = x1;
        j = x2;
    }
    x1path.reverse();
    x2path.reverse();
    (j,x1path,x2path)
}

#[inline]
pub fn reduce_by_half<D>(xs:&[D]) -> Vec<D> where D:Num+Copy{
    let two = D::one() + D::one();
    (0..xs.len()/2).map(|i|(xs[2*i]+xs[2*i+1])/two).collect()
}

fn expand_window(x1path:&[usize],x2path:&[usize],x2len:usize,x1len:usize,radius:usize)-> Vec<(usize,usize)>{
    debug_assert!(x1path.len() == x2path.len());
//    print_path(x1path,x2path);
    let mapped_path = x1path.iter().zip(x2path.iter()).map(|(&i,&j)|map_to_original(i,j))
        .fold(vec![],|mut acc,mut x|{acc.append(&mut x);acc});
    let mut region:Vec<Option<(usize,usize)>> = vec![None;x1len+1];
    for (i,j) in mapped_path{
        let start = if i > radius { i - radius } else { 1 };
        let end = if i + radius <= x1len { i + radius } else { x1len+1 };
        for k in start..end{
            let left = if j > radius { j - radius } else { 1 };
            let right = if j + radius <= x2len { j + radius }else{x2len+1};
            match region[k]{
                Some(prev) => region[k] = Some((prev.0,right)),
                None => region[k] = Some((left,right)),
            };
        }
    }
    (1..(x1len+1)).filter_map(|k|region[k].map(|range|(k,range)))
        .map(|(k,range)|get_window(k,range,x2len))
        .fold(vec![],|mut acc,mut x|{acc.append(&mut x);acc})
}

#[inline]
fn map_to_original(i:usize,j:usize)->Vec<(usize,usize)>{
    vec![(2*i,2*j),(2*i+1,2*j),(2*i,2*j+1),(2*i+1,2*j+1)]
}

#[inline]
fn get_window(k:usize,range:(usize,usize),x2len:usize)->Vec<(usize,usize)>{
    (range.0..range.1).filter(|&e| 0 < e && e <= x2len ).map(|e|(k,e)).collect()
}

// #[inline]
// fn print_window(window:&[(usize,usize)],x1len:usize,x2len:usize){
//     let x2start = window.iter().filter(|&&(i,_)| i <= 1 )
//         .map(|&(_,j)|j).min().unwrap();
//     print_window_sub(window,x1len,x2start);
// }

// #[inline]
// fn print_window_sub(window:&[(usize,usize)],x1len:usize,x2start:usize){
//     if x1len < 100 {
//         print_window_with_interval(window,x1len,x2start,1);
//     }else if x1len < 500 {
//         print_window_with_interval(window,x1len,x2start,5);
//     }else if x1len < 10000{
//         print_window_with_interval(window,x1len,x2start,10);
//     }else{
//         print_window_with_interval(window,x1len,x2start,100);
//     }
// }


// #[inline]
// fn print_window_with_interval(window:&[(usize,usize)],x1len:usize,x2start:usize,gap:usize){
//     use std::collections::HashSet;
//     let window:HashSet<_> = window.iter().clone().collect();
//     for i in 0..(x1len+1)/gap{
//             for j in x2start..x2start + (x1len+1)/gap{
//                 if window.contains(&(i*gap,j*gap)){
//                     eprint!("*");
//                 }else{
//                     eprint!("-");
//                 }
//             }
//         eprintln!();
//     }
// }



// #[inline]
// fn expand_line(i:usize,j:usize,x2len:usize,x1len:usize,radius:usize)->Vec<(usize,usize)>{
//     let mut res1 = expand_for_next_recursion(2*i,2*j,x1len,x2len,radius);
//     let mut res2 = expand_for_next_recursion(2*i+1,2*j+1,x1len,x2len,radius);
//     res1.append(&mut res2);
//     res1
// }

// #[inline]
// fn expand_for_next_recursion(i:usize,j:usize,x1len:usize,x2len:usize,radius:usize)->Vec<(usize,usize)>{
//     if i == 0 || i > x1len || j == 0 || j > x2len{
//         vec![]
//     }else{
//         let start = if j > radius + 1 { j - 1 - radius } else{ 1 };
//         let end = if j + radius <= x2len { j + radius } else { x2len + 1 };
//         (start .. end).map(|e|(i,e)).collect()
//     }
// }
