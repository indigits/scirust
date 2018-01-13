use std::cmp::{Ord, Ordering};

#[doc="Sort a list of keys and items with respect to the keys, in place using the quick sort algorithm.
"]
pub fn sort<TKey : Ord + Clone, TItem : Clone>(keys: &mut Vec<TKey>, items: &mut Vec<TItem>) {
    let count = keys.len();
    if count <= 1 {
        return;
    }

    if count == 2 {
        if keys[0].cmp(&keys[1]) == Ordering::Greater {
            swap(keys, 0, 1);
            swap(items, 0, 1);
        }
        return;
    }

    // insertion sort
    if count <= 10 {
        for i in 1..count {
            let key = keys[i].clone();
            let item = items[i].clone();
            let mut j = i - 1;
            while j >= 0 && keys[j].cmp(&key) == Ordering::Greater {
                keys[j + 1] = keys[j].clone();
                items[j + 1] = items[j].clone();
                j-=1;
            }
            keys[j + 1] = key;
            items[j + 1] = item;
        }
        return;
    }

    // local sort implementation
    quick_sort(keys, items, 0, count - 1);
}

#[doc="Recursive implementation for an in place quick sort on a list while reordering one other list accordingly.
keys - The list which is sorted using quick sort.
items - The list which is automatically reordered accordingly.
left - The left boundary of the quick sort.
right - The right boundary of the quick sort.
"]
fn quick_sort<T : Ord + Clone, TItems : Clone>(keys: &mut Vec<T>, items: &mut Vec<TItems>, left: usize, right: usize) {
    let mut left = left;
    let mut right = right;

    loop {
        // Pivoting
        let mut a = left;
        let mut b = right;
        let p = a + ((b - a) >> 1); // midpoint

        if keys[a].cmp(&keys[p]) == Ordering::Greater {
            swap(keys, a, p);
            swap(items, a, p);
        }

        if keys[a].cmp(&keys[b]) == Ordering::Greater {
            swap(keys, a, b);
            swap(items, a, b);
        }

        if keys[p].cmp(&keys[b]) == Ordering::Greater {
            swap(keys, p, b);
            swap(items, p, b);
        }

        let pivot = keys[p].clone();

        // Hoare Partitioning
        loop {
            while keys[a].cmp(&pivot) == Ordering::Greater {
                a+=1;
            }

            while pivot.cmp(&keys[b]) == Ordering::Greater {
                b-=1;
            }

            if a > b {
                break;
            }

            if a < b {
                swap(keys, a, b);
                swap(items, a, b);
            }

            a+=1;
            b-=1;

            if a > b {
                break;
            }
        }

        // In order to limit the recursion depth to log(n), we sort the
        // shorter partition recursively and the longer partition iteratively.
        if (b - left) <= (right - a) {
            if left < b {
                quick_sort(keys, items, left, b);
            }

            left = a;
        } else {
            if a < right {
                quick_sort(keys, items, a, right);
            }

            right = b;
        }

        if left >= right {
            break;
        }
    }
}

#[doc="Performs an in place swap of two elements in a list.
keys - The list in which the elements are stored.
a - The index of the first element of the swap.
b - The index of the second element of the swap.
"]
fn swap<T: Clone>(keys: &mut Vec<T>, a: usize, b: usize) {
    if a == b {
        return;
    }

    let local = keys[a].clone();
    keys[a] = keys[b].clone();
    keys[b] = local;
}