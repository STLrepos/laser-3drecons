#include <string>
#include <stdio.h>
#include <iostream>

#include "Laser3DReconst.hpp"

using namespace std;

int main(char* argc, char** argv)
{
  int ret =0 ;
  
  Laser3DReconst l3r;
  
  l3r.init();
  
  l3r.run();
  
  return ret;
}
