# building instructions

```sh
mkdir -p build && cd build            # create a build folder, if it does not exist
cmake .. -DCMAKE_BUILD_TYPE=Release   # configure cmake (only need to do this once), can also use =Debug (slower)
cmake --build .                       # build the program
cd ..                                 # exit the build folder
./build/grav-search                   # run the program
```
