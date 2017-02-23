#[macro_use]
extern crate bencher;
extern crate rand;

extern crate sralg;

// external dependencies
use rand::Rng;
use bencher::Bencher;

// local dependencies
use sralg::sort::*;
use sralg::sort::quick_sort::*;
use sralg::sort::insertion_sort::*;

fn bench_is_ascending_slice(b: &mut Bencher){
    let v = (0..10000).collect::<Vec<i32>>();
    b.iter(|| {
        assert!(is_ascending_slice(v.as_slice()));
            });
}
fn bench_is_ascending_buffer(b: &mut Bencher){
    let v = (0..10000).collect::<Vec<i32>>();
    b.iter(|| unsafe {
        assert!(is_ascending_buffer(v.as_ptr(), v.len()));
            });
}
fn bench_is_descending_slice(b: &mut Bencher){
    let v = (0..10000).map(|idx| (20000 - idx)).collect::<Vec<i32>>();
    b.iter(|| {
        assert!(is_descending_slice(v.as_slice()));
            });
}
fn bench_is_descending_buffer(b: &mut Bencher){
    let v = (0..10000).map(|idx| (20000 - idx)).collect::<Vec<i32>>();
    b.iter(|| unsafe {
        assert!(is_descending_buffer(v.as_ptr(), v.len()));
            });
}



fn qs_slice_reverse_data(b: &mut Bencher){
    let mut v = (0..1000).map(|idx| (20000 - idx)).collect::<Vec<i32>>();
    b.iter(|| {
        quick_sort_slice(v.as_mut_slice());
    });
}

fn qs_slice_random_data(b: &mut Bencher){
    // create a task-local Random Number Generator
    let mut rng = rand::thread_rng();
    let mut v: Vec<usize> = rng.gen_iter::<usize>().take(10000).collect();
    b.iter(|| {
        quick_sort_slice(v.as_mut_slice());
    });
}


fn qs_buffer_reverse_data(b: &mut Bencher){
    let mut v = (0..10000).map(|idx| (20000 - idx)).collect::<Vec<i32>>();
    b.iter(|| {
        unsafe {
            quick_sort_buffer(v.as_mut_ptr(), v.len());
        }
    });
}

fn qs_buffer_random_data(b: &mut Bencher){
    // create a task-local Random Number Generator
    let mut rng = rand::thread_rng();
    let mut v: Vec<usize> = rng.gen_iter::<usize>().take(10000).collect();
    b.iter(|| {
        unsafe {
            quick_sort_buffer(v.as_mut_ptr(), v.len());
        }
    });
}

benchmark_group!(qs, qs_slice_reverse_data, qs_slice_random_data, qs_buffer_reverse_data, qs_buffer_random_data);


fn ins_slice_reverse_data(b: &mut Bencher){
    let mut v = (0..10000).map(|idx| (20000 - idx)).collect::<Vec<i32>>();
    b.iter(|| {
        insertion_sort_slice(v.as_mut_slice());
    });
}

fn ins_slice_random_data(b: &mut Bencher){
    // create a task-local Random Number Generator
    let mut rng = rand::thread_rng();
    let mut v: Vec<usize> = rng.gen_iter::<usize>().take(10000).collect();
    b.iter(|| {
        insertion_sort_slice(v.as_mut_slice());
    });
}

fn ins_buffer_reverse_data(b: &mut Bencher){
    let mut v = (0..10000).map(|idx| (20000 - idx)).collect::<Vec<i32>>();
    b.iter(|| {
        unsafe {
            insertion_sort_buffer(v.as_mut_ptr(), v.len());
        }
    });
}

benchmark_group!(ins, ins_slice_reverse_data, ins_slice_random_data, ins_buffer_reverse_data);


benchmark_group!(benches, bench_is_ascending_slice, bench_is_ascending_buffer, bench_is_descending_slice, bench_is_descending_buffer);
benchmark_main!(benches, qs, ins);



