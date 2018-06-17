#include <fstream>
#include "CrystalLighting.h"

template <class T>
void getVector(vector<T>& v)
{
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

int main(int argc, char* argv[])
{
    CrystalLighting cl;
    int H;
    cin >> H;
    vector<string> targetBoard(H);
    getVector(targetBoard);
    int costLantern, costMirror, costObstacle, maxMirrors, maxObstacles;
    cin >> costLantern >> costMirror >> costObstacle >> maxMirrors >>
        maxObstacles;

    {
        std::ofstream of("/Users/tamas.kenez/tc/tco2018_2/lastgame");
        of << targetBoard.size() << "\n";
        for (auto& r : targetBoard) {
            of << r << "\n";
        }
        of << costLantern << "\n"
           << costMirror << "\n"
           << costObstacle << "\n"
           << maxMirrors << "\n"
           << maxObstacles << "\n";
    }

    vector<string> ret = cl.placeItems(targetBoard, costLantern, costMirror,
                                       costObstacle, maxMirrors, maxObstacles);
    cout << ret.size() << "\n";
    for (int i = 0; i < (int)ret.size(); ++i)
        cout << ret[i] << "\n";
    cout.flush();
    {
        cerr << ret.size() << "\n";
        for (int i = 0; i < (int)ret.size(); ++i)
            cerr << ret[i] << "\n";
        cerr.flush();
    }
}
