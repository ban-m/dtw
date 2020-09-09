use super::TRACE;
use num::Float;
/// dynamic time warping Itakura's pentagram
pub fn dtw_itakura<D, F, T>(
    _x1: &[D],
    _x2: &[D],
    _dist: &F,
    _band: usize,
) -> Result<(T, TRACE, usize), String>
where
    F: Fn(&D, &D) -> T,
    T: Float,
{
    // impliment someday.
    Ok((num::zero(), (vec![], vec![]), 0))
}
