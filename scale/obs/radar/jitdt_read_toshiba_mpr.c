/*
  JIT-DT version of TOSHIBA-formatted radar data reader
  Shigenori Otsuka
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <sys/time.h>
#include <jitclient.h>
#include "read_toshiba_mpr.h"
#include "jitdt_read_toshiba_mpr.h"

int jitdt_read_toshiba(int n_type, char *jitdt_place, mppawr_header hd[n_type],
                       float az[n_type][ELDIM][AZDIM], float el[n_type][ELDIM][AZDIM],
                       float rtdat[n_type][ELDIM][AZDIM][RDIM])
{
  const size_t bufsize = 40 * 1024 * 1024; // fixed size
  int i_type, ierr;
  int bsize[n_type];
  unsigned char *buf;
  char fname[(PATH_MAX + 1) * n_type - 1];
  struct timeval t0, t;
  const int opt_verbose=0;


  gettimeofday(&t0, NULL);

  buf = malloc(n_type * bufsize);
  if(buf == NULL){
    printf("failed to allocate memory in jitdt_read_toshiba");
    return -99;
  }

  for(i_type = 0; i_type < n_type; i_type++){
    bsize[i_type] = bufsize;
  }

  gettimeofday(&t, NULL);
  printf("......jitdt_read_toshiba:allocate_buffer:%15.6f\n", (float)(t.tv_sec-t0.tv_sec) + (float)(t.tv_usec-t0.tv_usec)/1000000.0);
  t0 = t;

  ierr = jitget(jitdt_place, fname, buf, bsize, n_type);

  gettimeofday(&t, NULL);
  printf("......jitdt_read_toshiba:jitget:%15.6f\n", (float)(t.tv_sec-t0.tv_sec) + (float)(t.tv_usec-t0.tv_usec)/1000000.0);
  t0 = t;

//  if(fname[0] == 0){
  if(ierr != 0) {
    printf("jitget failed: jitdt_place=%s\n", jitdt_place);
//    return -9;
    return ierr;
  }

  for(i_type = 0; i_type < n_type; i_type++){
    ierr = decode_toshiba_mpr(bsize[i_type], buf + i_type * bufsize, opt_verbose, hd + i_type, az[i_type], el[i_type], rtdat[i_type]);
    if(ierr != 0) return ierr;
  }

  gettimeofday(&t, NULL);
  printf("......jitdt_read_toshiba:decode_toshiba:%15.6f\n", (float)(t.tv_sec-t0.tv_sec) + (float)(t.tv_usec-t0.tv_usec)/1000000.0);
  t0 = t;

  free(buf);

  gettimeofday(&t, NULL);
  printf("......jitdt_read_toshiba:deallocate_buffer:%15.6f\n", (float)(t.tv_sec-t0.tv_sec) + (float)(t.tv_usec-t0.tv_usec)/1000000.0);
  t0 = t;

  return 0;
}
