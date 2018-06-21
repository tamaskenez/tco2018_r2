#include <array>
#include <cstdio>

using namespace std;
using A = array<uint8_t, 5>;
using B = array<A, 5>;

int main()
{
  B b;
  printf("%d\n", (int)sizeof(B));
  return 0;
}
