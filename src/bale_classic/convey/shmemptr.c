#include <shmem.h>
#include <stdio.h>
#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "mpp2shmem.h"

int
main(void)
{
  shmem_init();
  int* data = mpp_alloc(sizeof(int));
  shmem_barrier_all();

  if (shmem_my_pe() == 0) {
    void* friend = shmem_ptr(data, 1);
    *data = (friend == NULL);
    if (friend == NULL)
      fputs("ERROR: shmem_ptr() is not working.  Either fix your environment to\n"
            "make shmem_ptr() work, or configure --without-shmem-ptr.\n", stderr);
  }

  shmem_barrier_all();
  int status = shmem_int_g(data, 0);
  shmem_barrier_all();
  shmem_finalize();
  return status;
}
