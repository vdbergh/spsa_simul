#include <stdio.h>
#include "prng.h"

int main(){
  prng_t prng;
  prng_init(&prng);
  while(1){
    unsigned int prng32=prng>>32;
    size_t ret=fwrite(&prng32,4,1,stdout);
    if(ret==0){
      break;
    }
    prng_get(&prng);
  }
}
