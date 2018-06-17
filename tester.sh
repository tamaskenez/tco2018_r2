seed=$2
if [[ "$1" == "r" ]]; then
  cmake --build b --config Release
  java -cp dl CrystalLightingVis -exec "b/Release/fortester" -seed "$seed"
  # java -jar dl/tester.jar -exec "b/Release/fortester" -seed "$seed"
elif [[ "$1" == "d" ]]; then
  cmake --build b --config Debug
  java -cp dl CrystalLightingVis -exec "b/Debug/fortester" -seed "$seed"
  # java -jar dl/tester.jar -exec "b/Debug/fortester" -seed "$seed"
else
  echo "Missing \$1" <&2
  exit 1
fi
