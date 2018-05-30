/*
  JIT-DT version of TOSHIBA-formatted radar data reader
  Shigenori Otsuka
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <jitclient.h>
#include "read_toshiba.h"
#include "jitdt_read_toshiba.h"

int jitdt_read_toshiba(int n_type, char *jitdt_place, pawr_header hd[n_type],
                       float az[n_type][ELDIM][AZDIM], float el[n_type][ELDIM][AZDIM],
                       float rtdat[n_type][ELDIM][AZDIM][RDIM])
{
  const size_t bufsize = 40 * 1024 * 1024; // fixed size
  int i_type, ierr;
  int bsize[n_type];
  unsigned char *buf;
  char fname[(PATH_MAX + 1) * n_type - 1];

  buf = malloc(n_type * bufsize);
  if(buf == NULL){
    printf("failed to allocate memory in jitdt_read_toshiba");
    return -99;
  }

  for(i_type = 0; i_type < n_type; i_type++){
    bsize[i_type] = bufsize;
  }
  ierr = jitget(jitdt_place, fname, buf, bsize, n_type);

//  if(fname[0] == 0){
  if(ierr != 0) {
    printf("jitget failed: jitdt_place=%s\n", jitdt_place);
//    return -9;
    return ierr;
  }

  for(i_type = 0; i_type < n_type; i_type++){
    ierr = decode_toshiba(bsize[i_type], buf + i_type * bufsize, hd + i_type, az[i_type], el[i_type], rtdat[i_type]);
    if(ierr != 0) return ierr;
  }

  free(buf);

  return 0;
}
