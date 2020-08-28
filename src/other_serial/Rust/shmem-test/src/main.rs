use shmem::Shmem;
use shmem::ShmemAtomic;

fn main() {
    let s = Shmem::new().expect("Shmem initialization failed");
    let me = s.my_pe();
    let num = s.n_pes();
    //println!("Hello, world from PE {} of {}!", me, num);
    let mut v = s.new_obj::<i64>(100).expect("Allocation error");
    v.local_part[0] = 40;
    v.atomic_add(0,2,0).expect("atomic error");
    s.fence();
    if num > 1 && me == 0 {
        v.put(0,&[v.local_part[0]],2).expect("put error");
    }
    s.barrier();
    println!("At index 0 on {}/{} value is {}", me, num, v.local_part[0]);
}
