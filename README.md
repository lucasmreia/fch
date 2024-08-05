# A Fast High-Dimensional Continuation Hypercubes Algorithm

This repository contains the C++ implementation of the Fast Continuation Hypercubes algorithm.

This algorithm creates a piecewise linear approximation of a manifold of dimension $(n-k)$ which is implicitly defined by the function $F(\mathbf{x})=\mathbf{0}$, where $F:\mathbb{R}^n\to\mathbb{R}^k$.

For comparison, we also provide implementations for the permutahedron-based tracing algorithm of Boissonnat et al. [1], and the Generalized Combinatorial Continuation Hypercubes [2].

https://github.com/user-attachments/assets/ab187621-dc49-452d-ad52-46267e2ed2b4

## License

The files in this repository are free to use as long as you cite, or give credit to, the original creators.

## Requirements
 - A C++ compiler with C++20 support (tested on GCC and Apple CLang).
 - cmake
 - git

## Dependencies
 - [Eigen](https://gitlab.com/libeigen/eigen) 
 - [argparse](https://github.com/p-ranav/argparse)
 - [zp7](https://github.com/zwegner/zp7) (used when PDEP and PEXT are missing)

## Implemented algorithms
 - Permutahedron-based Tracing Algorithm (PTA)
 - Generalized Combinatorial Continuation Hypercubes (GCCH)
 - Fast Continuation Hypercubes (FCH)

## Building and Running

### 1. Cloning the repository

  ```sh
  git clone --recursive https://github.com/lucasmreia/fch
  cd fch
  ```

### 2. Compiling

The test case is selected by a definition in the build command. For example, to compile a program for the Klein Bottle test case ($F:\mathbb{R}^5\to\mathbb{R}^3$), you use:

  ```
  -DTEST_CASE=KLEIN_5TO3
  ```

The available test cases are:

  ```
  enum class TestCase {
      KLEIN_5TO3,
      HYPERSHPERE_5TO1,
      CIRCLE_9TO8,
  };
  ```

#### Preparing the build

That being said, create the directory for the build files and prepare the build:

  ```sh
  mkdir build
  cmake -DCMAKE_BUILD_TYPE=Release -DTEST_CASE=KLEIN_5TO3 -S . -B build
  ```

If you want to change the test case, you can rerun the second command with a different case:

  ```sh
  cmake --build build --target clean
  cmake -DCMAKE_BUILD_TYPE=Release -DTEST_CASE=HYPERSHPERE_5TO1 -S . -B build
  ```

#### Compiling the PTA

  ```sh
  cmake --build build --target pta skeleton_pta -- -j 4
  ```

#### Compiling the GCCH

  ```sh
  cmake --build build --target gcch reorder_hypercubes skeleton_hypercubes -- -j 4
  ```

#### Compiling the FCH

  ```sh
  cmake --build build --target fch fch_to_hypercubes reorder_hypercubes skeleton_hypercubes -- -j 4
  ```


### 3. Running

Create the directory for the output files:

  ```sh
  mkdir out
  ```

The usual pipeline is first to generate the binary output of an algorithm, then convert it to the .pol format for visualization. Keep in mind that if any of the algorithms do not produce any output, it means that the coordinates of the first point (the starting point for the continuation method) must be tweaked.

The Combinatorial Skeleton[2] is used to create the high-dimensional cells of the output, producing the .pol file.

The pipeline was split into multiple executables for convenience. You can check the help of each executable individually.

#### Running the PTA

  ```sh
  ./build/pta -o out/pta.bin
  ./build/skeleton_pta -i out/pta.bin -o out/pta.pol -ff 15.6
  ```

#### Running the GCCH

  ```sh
  ./build/gcch -o out/gcch.bin
  ./build/reorder_hypercubes -i out/gcch.bin -o out/gcch_reordered.bin
  ./build/skeleton_hypercubes -i out/gcch_reordered.bin -o out/gcch.pol -ff 15.6
  ```

#### Running the FCH

  ```sh
  ./build/fch -o out/fch.bin
  ./build/fch_to_hypercubes -i out/fch.bin -o out/fch_hypercubes.bin
  ./build/reorder_hypercubes -i out/fch_hypercubes.bin -o out/fch_hypercubes_reordered.bin
  ./build/skeleton_hypercubes -i out/fch_hypercubes_reordered.bin -o out/fch.pol -ff 15.6
  ```

## References

[1] Boissonnat, J. D., Kachanovich, S., & Wintraecken, M. (2023). Tracing Isomanifolds in d in Time Polynomial in d using Coxeter–Freudenthal–Kuhn Triangulations. SIAM Journal on Computing, 52(2), 452-486.

[2] Castelo, A., Nakassima, G., Bueno, L.M. et al. A generalized combinatorial marching hypercube algorithm. Comp. Appl. Math. 43, 127 (2024). https://doi.org/10.1007/s40314-024-02627-4
