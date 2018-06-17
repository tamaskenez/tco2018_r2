#include <fstream>
#include "CrystalLighting.h"

int main()
{
    std::ifstream f("/Users/tamas.kenez/tc/tco2018_2/lastgame");
  assert(f.good());
    int NR;
    f >> NR;
    vector<string> tb;
    tb.reserve(NR);
    for (int i = 0; i < NR; ++i) {
        string r;
        f >> r;
        tb.emplace_back(move(r));
    }
    int CL, CM, CO, MM, MO;
    f >> CL >> CM >> CO >> MM >> MO;
    CrystalLighting cl;
    vector<string> ret = cl.placeItems(tb, CL, CM, CO, MM, MO);
    return EXIT_SUCCESS;
}
