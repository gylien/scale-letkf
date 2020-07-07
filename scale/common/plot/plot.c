#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image.h"

int plot(int nwidth,int nheight, char *idata, int nchar, char *cfile) {
  int i;
  image_t *img = NULL;
  uint8_t x,y;
  char cfilename[200]; 
  int imode = 2;

  if ( nwidth > 600 ) {
     printf("plot.c : abort - nwidth should be smaller than 600");
	  return 99;
}

  if ( (img = allocate_image(nwidth, nheight, COLOR_TYPE_INDEX)) == NULL ) {
    printf("plot.c : abort error in allocate_image");          
	  return 99;
  }
 
/* rain rate color scale for weather.riken.jp */  
  img->palette_num = 12;
  img->palette[0]  = color_from_rgb(128, 128, 128);
  img->palette[1]  = color_from_rgb(255, 255, 255);
  img->palette[2]  = color_from_rgb(200, 171, 216);
  img->palette[3]  = color_from_rgb(144,  86, 177);
  img->palette[4]  = color_from_rgb( 57,   0, 177);
  img->palette[5]  = color_from_rgb(  0,   7, 246);
  img->palette[6]  = color_from_rgb(  5, 125,  87);
  img->palette[7]  = color_from_rgb( 98, 207,   9);
  img->palette[8]  = color_from_rgb(255, 253,   0);
  img->palette[9]  = color_from_rgb(255, 190,   0);
  img->palette[10] = color_from_rgb(255,  98,   0);
  img->palette[11] = color_from_rgb(255,   0,   0);

  for (y = 0; y < img->height; y++) {
    for (x = 0; x < img->width; x++) {
      pixcel_t *p = &img->map[y][x];
      p->i = (unsigned char)idata[y*nwidth+x];
    }
  }
  
  img=image_to_rgba(img);
  for (y = 0; y < img->height; y++) {
    for (x = 0; x < img->width; x++) {
      pixcel_t *p = &img->map[y][x];
      if (idata[y*nwidth+x] == 1){
        p->c.a = 0;
      }else{
        p->c.a = 216;
      }
    }
  }

  strncpy(cfilename, cfile, nchar);
  cfilename[nchar]='\0';
 
  if (img != NULL) {
      write_png_file(cfilename, img);
      free_image(img);
      img = NULL;
    } 
  return 0;
}
