use super::*;
#[test]
fn test_of_test(){
    assert!(true);
}

#[test]
fn test_dtw(){
    let x1 = vec![1.,1.];
    let x2 = vec![2.,2.];
    let res = dtw::normal_dtw(&x1,&x2);
    debug_assert!((res-2.0).abs() < 0.001,"{}",dtw::normal_dtw(&x1,&x2));
}
#[test]
fn test_dtw2(){
    let x1 = vec![1.];
    let x2 = vec![1.];
    debug_assert!(dtw::normal_dtw(&x1,&x2) < 0.001,"{}",dtw::normal_dtw(&x1,&x2));
}

#[test]
fn test_dtw3(){
    let x = vec![1.,2.,3.,4.,5.];
    let y = vec![2.,3.,4.];
    let res = dtw::normal_dtw(&x,&y);
    debug_assert!((res - 2.).abs() < 0.1,"{}",res)
}
#[test]
fn test_dtw4(){
    let x = vec![1., 1., 1., 2., 2., 2., 2.];
    let y = vec![1., 1., 2.,2.];
    let res = dtw::subdtw(&y,&x);
    println!("{}",res);
    debug_assert!(res < 0.001,"{}",res);
}
#[test]
fn test_dtw5(){
    let x = vec![1., 1., 1.,1.,1.];//,2., 2., 2., 2., 3., 2., 0.];
    let y = vec![2.,0., 0., 1.];//, 1., 2., 4., 2., 1., 2., 0.];
    let res = dtw::subdtw(&x,&y);
    debug_assert!(res -3.0 < 0.001,"{}",res);
}

#[test]
fn test_dtw6(){
    let x = vec![1., 1., 4.,5.,1.];//,2., 2., 2., 2., 3., 2., 0.];
    let y = vec![6.,6.];//, 1., 2., 4., 2., 1., 2., 0.];
    let res = dtw::subdtw(&y,&x);
    debug_assert!(res -3.0 < 0.001,"{}",res);
    debug_assert!(res > 0.,"{}",res);
}

#[test]
fn test_dtw7(){
    let x = vec![1., 1., 4.,5.,1.];
    let y = vec![6.,6.];
    let d = |x:&f32,y:&f32| (x-y).powi(2);
    let (res,(ypath,xpath),start) = dtw::dtw(&y,&x,dtw::Mode::Sub,&d).unwrap();
    println!("{},{}\n{:?}\n{:?}",start,res,xpath,ypath);
    debug_assert!(start == 3 ,"{},{}\n{:?}\n{:?}",start,res,xpath,ypath);
}

#[test]
fn test_dtw8(){
    let _ = normalize::NormalizeType::Z;
    let x = vec![1., 1., 4.,5.,1.];//,2., 2., 2., 2., 3., 2., 0.];
    let y = vec![6.,6.];//, 1., 2., 4., 2., 1., 2., 0.];
    let (res,(ypath,xpath)) = dtw::subdtw_with_path(&y,&x);
    assert_eq!(xpath.len(),ypath.len());
    println!("{}",res);
    println!("{:?}\n{:?}",x,y);
    for yi in ypath{
        print!("{:3} ",y[yi]);
    }
    println!("");
    for xi in xpath{
        print!("{:3} ",x[xi]);
    }
    println!("");
    let _ = normalize::NormalizeType::Z;
    let x = vec![1., 1., 3.,1.,1.];//,2., 2., 2., 2., 3., 2., 0.];
    let y = vec![2.,0., 3., 1.];//, 1., 2., 4., 2., 1., 2., 0.];
    let (res,(ypath,xpath)) = dtw::subdtw_with_path(&y,&x);
    assert_eq!(xpath.len(),ypath.len());
    println!("{}",res);
    println!("{:?}\n{:?}",x,y);
    for yi in ypath{
        print!("{:3} ",y[yi]);
    }
    println!("");
    for xi in xpath{
        print!("{:3} ",x[xi]);
    }
    println!("");
    assert!(true);
}

#[test]
fn test_fastdtw_small(){
    let radius = 16;
    let num = 16;
    let x:Vec<_> = (0..num).map(|i:usize| (2. * i as f32 * std::f32::consts::PI/ num as f32)).collect();
    let y:Vec<_>= (0..num).map(|i:usize| (2. * i as f32 * std::f32::consts::PI/ num as f32)).collect();
    let dist = |x:&f32,y:&f32| (x-y).powi(2);
    let (score,_,_) = fast_dtw(&x,&y,&dist,radius,false).unwrap();
    eprintln!("{}",score);
    debug_assert!(score<0.01,"{}",score);
}
#[test]
fn test_fastdtw_long(){
    let radius = 16;
    let num = 256;
    let x:Vec<_> = (0..num).map(|i:usize| (2. * i as f32 * std::f32::consts::PI/ num as f32)).collect();
    let y:Vec<_>= (0..num).map(|i:usize| (2. * i as f32 * std::f32::consts::PI/ num as f32)).collect();
    let dist = |x:&f32,y:&f32| (x-y).powi(2);
    let (score,_,_) = fast_dtw(&x,&y,&dist,radius,false).unwrap();
    eprintln!("{}",score);
    debug_assert!(score<0.01,"{}",score);
}

#[test]
fn test_fastdtw_long_not_powerof2(){
    let radius = 30;
    let num = 600;
    let x:Vec<_> = (0..num).map(|i:usize| (2. * i as f32 * std::f32::consts::PI/ num as f32)).collect();
    let y:Vec<_>= (0..num).map(|i:usize| (2. * i as f32 * std::f32::consts::PI/ num as f32)).collect();
    let dist = |x:&f32,y:&f32| (x-y).powi(2);
    let (score,_,_) = fast_dtw(&x,&y,&dist,radius,false).unwrap();
    eprintln!("{}",score);
    debug_assert!(score<0.01,"{}",score);
}

#[test]
fn test_fastdtw_sub_small(){
    let radius = 16;
    let num = 16;
    let y:Vec<_> = (0..3 * num).map(|i:usize|if num < i && i < 2*num {
        (2. * (i-num) as f32 * std::f32::consts::PI/ num as f32)
    }else{
        0.
    }).collect();
    let x:Vec<_>= (0..num).map(|i:usize| (2. * i as f32 * std::f32::consts::PI/ num as f32)).collect();
    let dist = |x:&f32,y:&f32| (x-y).powi(2);
    let (score,(x1path,x2path),_) = fast_dtw(&x,&y,&dist,radius,true).unwrap();
    for (&i,&j) in x1path.iter().zip(x2path.iter()){
        eprint!("({},{})->",i,j);
    }
    eprintln!("finish");
    debug_assert!(score<0.01,"{}",score);

}

#[test]
fn test_fastdtw_sub_large(){
    let radius = 16;
    let num = 32;
    let y:Vec<_> = (0..3 * num).map(|i:usize|if num < i && i < 2*num {
        (2.*(i-num) as f32 * std::f32::consts::PI/ num as f32).sin()
    }else{
        0.
    }).collect();
    let x:Vec<_>= (0..num).map(|i:usize| (2. * i as f32 * std::f32::consts::PI/ num as f32).sin()).collect();
    eprintln!("{:?}",y);
    let y1 = fastdtw::reduce_by_half(&y);
    let x1 = fastdtw::reduce_by_half(&x);
    eprintln!("{:?}",y1);
    let dist = |x:&f32,y:&f32| (x-y).powi(2);
    let (score,(x1path,x2path),_) = fast_dtw(&x1,&y1,&dist,radius,true).unwrap();
    for (&i,&j) in x1path.iter().zip(x2path.iter()){
        eprint!("({},{})->",i,j);
    }
    eprintln!("finish");
    debug_assert!(score<0.01,"{}",score);
}

#[test]
fn test_fastdtw_sub_large_not_pow2(){
    let radius = 20;
    let num = 400;
    let y:Vec<_> = (0..3 * num).map(|i:usize|if num < i && i < 2*num {
        (2.*(i-num) as f32 * std::f32::consts::PI/ num as f32).sin()
    }else{
        0.
    }).collect();
    let x:Vec<_>= (0..num).map(|i:usize| (2. * i as f32 * std::f32::consts::PI/ num as f32).sin()).collect();
    eprintln!("{:?}",y);
    let y1 = fastdtw::reduce_by_half(&y);
    let x1 = fastdtw::reduce_by_half(&x);
    eprintln!("{:?}",y1);
    let dist = |x:&f32,y:&f32| (x-y).powi(2);
    let (_,(x1path,x2path),_) = fast_dtw(&x1,&y1,&dist,radius,true).unwrap();
    for (&i,&j) in x1path.iter().zip(x2path.iter()){
        eprint!("({},{})->",i,j);
    }
    eprintln!("finish");
    let (score,_,_) = fast_dtw(&x,&y,&dist,radius,true).unwrap();
    debug_assert!(score<0.01,"{}",score);
}
#[test]
fn test_how_good(){
    use rand::thread_rng;
    use rand::Rng;
    let mut rng = thread_rng();
    let radius = 20;
    let num = 20 ;
    let dist = |x:&f32,y:&f32| (x-y).powi(2);
    let x :Vec<f32> =  rng.gen_iter().take(num).collect();
    let y :Vec<f32> = rng.gen_iter().take(num).collect();
    let (score,_,_) = fast_dtw(&x,&y,&dist,radius,true).unwrap();
    assert!(score != 0.);
}
#[test]
fn test_how_good2(){
    use rand::thread_rng;
    use rand::Rng;
    let mut rng = thread_rng();
    let radius = 200;
    let num = 2000;
    let dist = |x:&f32,y:&f32| (x-y).powi(2);
    let x :Vec<f32> =  rng.gen_iter().take(num).collect();
    let y :Vec<f32> = rng.gen_iter().take(num).collect();
    let (score,_,_) = fast_dtw(&x,&y,&dist,radius,false).unwrap();
    let (optscore,_,_) = dtw::dtw(&x,&y,dtw::Mode::Full,&dist).unwrap();
    debug_assert!((score-optscore)/optscore < 0.05,"{},{}",score,optscore);
}

// fn test_how_good_ab(reference:usize,query:usize,radius:usize){
//     use rand::thread_rng;
//     use rand::Rng;
//     let mut rng = thread_rng();
//     let dist = |x:&f32,y:&f32| (x-y).powi(2);
//     let x :Vec<f32> =  rng.gen_iter().take(query).collect();
//     let y :Vec<f32> = rng.gen_iter().take(reference).collect();
//     let (score,_,_) = fast_dtw(&x,&y,&dist,radius,true).unwrap();
//     let (optscore,_,_) = dtw::dtw(&x,&y,dtw::Mode::Sub,&dist).unwrap();
//     eprintln!("(r:{},q:{}r:{}){}",reference,query,radius,(score-optscore)/optscore);
//     assert!(true);
// }

// #[test]
// fn test_how_good3(){
//     let queries = vec![500];
//     let refsizes:Vec<_> = vec![1000].into_iter().map(|e|e*1000).collect();
//     let radius = vec![10,20,30,40,50,100,200,400];
//     for q in queries{
//         for re in refsizes.iter(){
//             for r in radius.iter(){
//                 test_how_good_ab(*re,q,*r);
//             }
//         }
//     }
//     assert!(false);
// }


#[test]
fn quick(){
    let mut reference = noise(100);
    reference.append(&mut noisy_courve(100));
    reference.append(&mut noise(2000));
    for _ in 0..10{
        let query = noisy_courve(200);
        if let Ok((res,_,_)) = super::dtw::dtw(&query,&reference,super::dtw::Mode::Sub,&hill){
            if let Ok((res2,_,_)) = super::dtw::dtw(&query,&reference,super::dtw::Mode::QuickSub,&hill){
                assert!((res-res2)<0.1,"{},{}",res,res2)
            }else{
                assert!(false)
            }
        }else{
            assert!(false)
        }
    }
}


fn sin_curve(t:usize)->Vec<f32>{
    (0..t).map(|i| (((2*i) as f32)/t as f32 * std::f32::consts::PI).sin()).collect()
}
fn noise(t:usize)->Vec<f32>{
    use rand::thread_rng;
    use rand::distributions::{IndependentSample,Normal};
    let mut rng = thread_rng();
    let normal = Normal::new(0.,1.,);
    (0..t).map(|_| normal.ind_sample(&mut rng) as f32).collect()
}

fn get_score(querysize:usize,refsize:usize,
             radius:usize,n:i32,alpha:f32)->Option<f64>{
    let hill = move |x:&f32,y:&f32|{
        let d = (x-y).powi(2*n);
        d/(alpha + d)
    };
    let query = sin_curve(querysize);
    let reference = noise(refsize);
    fast_dtw(&query,&reference,&hill,radius,true)
        .map(|e|e.0 as f64).ok()
}

#[test]
fn test_dp(){
    eprintln!();
    get_score(50,100,2,1,0.1);
}

#[test]
fn manytimes(){
    let querysizes = 500;
    let refsize = 100_000;
    for _ in 0..1{
        if let Some(res) = get_score(querysizes,refsize,50,3,0.74){
            assert!(!res.is_infinite());
        }
    }
}

fn read_eve(path:&std::path::Path)->Result<(Vec<f32>,String),()>{
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;
    let pathname = match path.to_str().map(|e|e.to_string()){
        Some(res) => res,
        None => return Err(()),
    };
    let f = File::open(path).map_err(|_|())?;
    let query = &BufReader::new(f)
        .lines().skip(1).filter_map(|e|e.ok())
        .filter_map(|line|line.split(',').nth(3)
                    .and_then(|e| e.parse().ok()))
        .take(1000)
        .collect();
    Ok((normalize(&query,NormalizeType::Z),pathname))
}

#[test]
fn signal_test(){
    use std::fs::DirEntry;
    use std::fs::read_dir;
    use std::path::Path;
    let hill =  |x:&f32,y:&f32|{
        let d = (x-y).powi(6);
        d/(0.01 + d)
    };
    let queries :Vec<_>= read_dir(&Path::new("./src/testdata/")).unwrap()
        .filter_map(|e|e.ok())
        .map(|entry:DirEntry| entry.path())
        .filter_map(|path|read_eve(&path).ok())
        .collect();
    for &(ref q1,_) in queries.iter(){
        for &(ref q2,_) in queries.iter(){
            let fast = dtw::dtw(q1,q2,dtw::Mode::FastSub(50),&hill).unwrap();
            let full = dtw::dtw(q1,q2,dtw::Mode::Sub,&hill).unwrap();
            assert!(fast.0 == 0. || (fast.0-full.0)/full.0*100. < 10.
                    ,"{}\t{}\t{}",fast.0,full.0,(fast.0-full.0)/full.0*100.);
        }
    }
}

#[inline]
fn hill(x:&f32,y:&f32)->f32{
    let d = (x-y).powi(8);
    d/(1.+d)
}

#[test]
fn mock(){
    use rand::distributions::{Normal,IndependentSample};
    use rand::thread_rng;
    let mut rng = thread_rng();
    let radius = 50;
    let normal = Normal::new(0.,1.,);
    let y:Vec<_> = (0..1000).map(|_| normal.ind_sample(&mut rng) as f32).collect();
    for _ in 0..9{
        let x:Vec<_> = (0..500).map(|_| normal.ind_sample(&mut rng) as f32).collect();
        let (score,_,_) = fast_dtw(&x,&y,&hill,radius,true).unwrap();
        let (optscore,_,_) = dtw::dtw(&x,&y,dtw::Mode::Sub,&hill).unwrap();
        eprintln!("{}\t{}\t{}",score,optscore,(score-optscore)/optscore);
    }
    assert!(true);
}

fn noisy_courve(t:usize)->Vec<f32>{
    use rand::distributions::{Normal,IndependentSample};
    use rand::thread_rng;
    let mut rng = thread_rng();
    let normal = Normal::new(0.,0.5,);
    sin_curve(t).into_iter()
        .map(|e| e + normal.ind_sample(&mut rng) as f32).collect()
}

#[test]
fn noisy(){
    use rand::distributions::{Normal,IndependentSample};
    use rand::thread_rng;
    let mut rng = thread_rng();
    let mut rng2 = thread_rng();
    let radius = 50;
    let normal = Normal::new(0.,0.5,);
    eprintln!();
    for _ in 0..1{
        let y:Vec<_> = (0..450).map(|_|normal.ind_sample(&mut rng) as f32)
            .chain(noisy_courve(100).into_iter())
            .chain((0..450).map(|_|normal.ind_sample(&mut rng2)as f32))
            .collect();
        let x:Vec<_> = noisy_courve(180);
        let (score,_,l1) = fast_dtw(&x,&y,&hill,radius,true).unwrap();
        let (optscore,_,l2) = dtw::dtw(&x,&y,dtw::Mode::Sub,&hill).unwrap();
        for p in y {
            eprintln!("{}",p);
        }
        eprintln!();
        for p in x{
            eprintln!("{}",p);
        }
        // assert!((score-optscore)/optscore<0.05,"{}\t{}\t{}\t:{}\t{}",score,optscore,(score-optscore)/optscore,l1,l2);
        eprintln!("{}\t{}\t{}\t:{}\t{}",score,optscore,(score-optscore)/optscore,l1,l2);
    }
}

#[test]
fn should_false(){
    let radius = 50;
    let querylen = 1000;
    let referencelen = 10_000;
    let mut short_noise = noise(querylen);
    let query :Vec<_>= sin_curve(querylen).iter()
        .zip(short_noise.iter())
        .map(|(&s,&n)| s + n)
        .collect();
    short_noise.append(&mut noise(referencelen-querylen));
    let mut reference = vec![0.;referencelen/2-querylen/2];
    reference.append(&mut sin_curve(querylen));
    reference.append(&mut vec![0.;referencelen/2-querylen/2]);
    let (score,_,l1) = fast_dtw(&query,&reference,&hill,radius,true).unwrap();
    let (optscore,_,l2) = dtw::dtw(&query,&reference,dtw::Mode::Sub,&hill).unwrap();
    assert!(true,"{} from {} vs {} from {}",score,l1,optscore,l2);
}

#[test]
fn test_how_good6(){
    use rand::thread_rng;
    use rand::Rng;
    let mut rng = thread_rng();
    let radius = 50;
    let reference = 1_00;
    let query = 500;
    let dist = |x:&f32,y:&f32| (x-y).powi(2);
    let x :Vec<f32> = sin_curve(query);
    let xs:Vec<f32> = sin_curve(query*10/18);
    let mut noise:Vec<f32> = sin_curve(reference);
    rng.shuffle(&mut noise);
    let y :Vec<f32> = noise.iter().chain(xs.iter()).map(|&e|e).collect();
    rng.shuffle(&mut noise);
    let y :Vec<f32>= y.into_iter().chain(noise.into_iter()).collect();
    let (score,(x1,y1),_) = fast_dtw(&x,&y,&dist,radius,true).unwrap();
    let (optscore,(x2,y2),_) = dtw::dtw(&x,&y,dtw::Mode::Sub,&dist).unwrap();
    eprintln!("fast:{}",score);
    for (i,j) in x1.into_iter().zip(y1.into_iter()){
        eprint!("({},{})->",i,j);
    }
    eprintln!("finish");
    eprintln!("sub:{}",optscore);
    for (i,j) in x2.into_iter().zip(y2.into_iter()){
        eprint!("({},{})->",i,j);
    }
    eprintln!("finish");
    assert!((score-optscore)/optscore < 0.05,"{},{}",score,optscore);
}

#[test]
fn test_ucr(){
    let mut reference = noise(1000);
    reference.append(&mut noisy_courve(200));
    reference.append(&mut noise(10000));
    let bandwidth = 20;
    let _ = ucr_dtw::DynamicTimeWarping::new(reference.clone(),hill,bandwidth);
    assert!(true);    
}
#[test]
fn mini_ucr(){
    let reference = noise(100);
    let bandwidth = 49;
    let dtw = ucr_dtw::DynamicTimeWarping::new(reference.clone(),hill,bandwidth);
    let query = noise(100);
    let (tes,_) = dtw.dtw(&query);
    let (score,_,_) = super::dtw::dtw(&query,&reference,super::dtw::Mode::Full,&hill).unwrap();
    assert!((score-tes).abs() < 0.01,"{},{}",score,tes);
}

#[test]
fn sin_ucr(){
    let reference = sin_curve(100);
    let bandwidth = 9;
    let dtw = ucr_dtw::DynamicTimeWarping::new(reference.clone(),hill,bandwidth);
    let query = sin_curve(100);
    let (tes,_) = dtw.dtw(&query);
    let (score,_,_) = super::dtw::dtw(&query,&reference,super::dtw::Mode::Sub,&hill).unwrap();
    assert!((score-tes).abs() < 0.01,"{},{}",score,tes);
}

#[test]
fn chiba_ucr(){
    let reference = sin_curve(100);
    let bandwidth = 9;
    let dtw = ucr_dtw::DynamicTimeWarping::new(reference.clone(),hill,bandwidth);
    let query = sin_curve(100);
    let (tes,_) = dtw.dtw(&query);
    let (score,_,_) = super::dtw::dtw(&query,&reference,super::dtw::Mode::SakoeChiba(9),&hill).unwrap();
    assert!((score-tes).abs() < 0.01,"{},{}",score,tes);
}
#[test]
fn exec_ucr(){
    let mut reference = noise(1000);
    reference.append(&mut sin_curve(200));
    reference.append(&mut noise(10000));
    let bandwidth = 50;
    let dtw = ucr_dtw::DynamicTimeWarping::new(reference.clone(),hill,bandwidth);
    let query = sin_curve(200);
    let (score,pos) = dtw.dtw(&query);
    let (tes,_,position) = super::dtw::dtw(&query,&reference,super::dtw::Mode::Sub,&hill).unwrap();
    assert!((score-tes).abs() < 0.01,"{},{},{},{}",score,tes,pos,position);
}


#[test]
fn threshold(){
    let reference = noise(1000);
    for _ in 0..10 {
        let query = noise(100);
        let score = dtw::dtw(&query,&reference,dtw::Mode::QuickSub,&hill).unwrap().0;
        let score2 = dtw::dtw(&query,&reference,dtw::Mode::Sub,&hill).unwrap().0;
        let thre = thresholddtw(&query,&reference,&hill,100000.).unwrap().0;
        assert!((score-thre).abs()<0.01,"{},{},{}",score,thre,score2);
    }
}

#[test]
fn threshold_small(){
    let reference = vec![1., 1., 1., 2., 2., 2., 2.];
    let query = vec![1., 1., 2.,2.];
    let score = dtw::dtw(&query,&reference,dtw::Mode::Sub,&hill).unwrap().0;
    let thre = thresholddtw(&query,&reference,&hill,100000.).unwrap().0;
    assert!((score-thre).abs()<0.01,"{},{}",score,thre);
}

#[test]
fn scouting_small(){
    let query = noisy_courve(500);
    let reference:Vec<_> = noise(500).into_iter()
        .chain(noisy_courve(500).into_iter())
        .chain(noise(1000).into_iter())
        .collect();
    let score = dtw::dtw(&query,&reference,dtw::Mode::Sub,&hill).unwrap();
    let scout = scouting_dtw(&query,&reference,&hill,Some(3),Some(1)).unwrap();
    assert!((score.0-scout.0).abs() < 0.01,"{},{}",score.0,scout.0);
}

#[test]
fn scouting_manytimes(){
    let times = 100;
    let sum :f32= (0..times).map(|_|{
        let query = noisy_courve(500);
        let reference:Vec<_> = noise(500).into_iter()
            .chain(noisy_courve(500).into_iter())
            .chain(noise(1000).into_iter())
            .collect();
        let score = dtw::dtw(&query,&reference,dtw::Mode::Sub,&hill).unwrap();
        let scout = scouting_dtw(&query,&reference,&hill,None,None).unwrap();
        (score.0 -scout.0).abs()/score.0})
        .sum();
    // assertion pass when overall error rate is less than 0.2 percent.
    assert!(sum <= 0.2,"average error:{}",sum *100./ times as f32);
}
#[test]
fn scouting_when_maxpack(){
    let reference:Vec<_> = noise(500);
    let query = noise(500);
    let sub = dtw::dtw(&query,&reference,dtw::Mode::Sub,&hill).unwrap().0;
    let scout = scouting_dtw(&query,&reference,&hill,Some(1),Some(1)).unwrap().0;
    assert!(sub == scout,"{},{}",sub,scout);
}

#[test]
fn scouting_maxpack_manytimes(){
    let times = 100;
    let sum :f32= (0..times).map(|_|{
        let query = noisy_courve(600);
        let reference:Vec<_> = noise(500).into_iter()
            .chain(noisy_courve(500).into_iter())
            .chain(noise(1000).into_iter())
            .collect();
        let score = dtw::dtw(&query,&reference,dtw::Mode::Sub,&hill).unwrap();
        let scout = scouting_dtw(&query,&reference,&hill,Some(1),Some(1)).unwrap();
        (score.0 -scout.0).abs()/score.0})
        .sum();
    assert!(sum <= 0.01,"average error:{}",sum *100./ times as f32);
}
