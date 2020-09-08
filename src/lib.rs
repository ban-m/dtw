#![crate_name = "dtw"]
#![crate_type = "lib"]
#![warn(missing_docs)]

//!A tiny implementation of dynamic time warping for Rust language.
//!
extern crate num;
extern crate rand;
extern crate order_stat;
mod dtw;
mod dtw_chiba;
mod dtw_itakura;
mod fastdtw;
mod result;
mod normalize;
mod quickdtw;
mod ucr_dtw;
mod thresholddtw;
mod scoutingdtw;
mod nw;
/// module for utility such as optimal dynamic time warping.
/// also some other convinient functions are here.
pub mod utils;
pub use dtw::*;
pub use fastdtw::fast_dtw;
pub use normalize::NormalizeType;
pub use normalize::normalize;
pub use normalize::normalize_mut;
pub use normalize::histgram_equalization;
pub use normalize::histgram_modify;
pub use ucr_dtw::DynamicTimeWarping;
pub use thresholddtw::thresholddtw;
pub use scoutingdtw::{scouting_dtw,scouting_threshold_dtw};
pub use nw::needleman_dtw;
#[cfg(test)]
mod test;
