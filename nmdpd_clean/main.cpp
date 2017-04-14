#include "NMDPD.h"
using namespace NMDPD_NS;

int main(int argc, char** argv)
{
  NMDPD* nm = new NMDPD(argc, argv);
  nm->exec();
  return 0;
}
