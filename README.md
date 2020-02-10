# DCEGM


[![Build Status](https://travis-ci.org/floswald/DCEGM.jl.svg?branch=master)](https://travis-ci.org/floswald/DCEGM.jl)

[![Build status](https://ci.appveyor.com/api/projects/status/dxcqu2mfiskgw90m?svg=true)](https://ci.appveyor.com/project/floswald/dcegm-jl)


Julia implementation of the matlab version at [https://github.com/fediskhakov/dcegm](https://github.com/fediskhakov/dcegm).

## Performance: `x7.5`

```julia
julia> DCEGM.bm()
Hi! This is Matlab version 9.7.0.1261785 (R2019b) Update 3 running on my laptop
t: 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
Retirement model solved with
500 asset points
25 periods at
0.0000020 lambda  
0.350 sigma
in  1.509s
wrote policy and value function to ascii in output/ . exiting matlab.
julia timing:
    0.206208 seconds (220.33 k allocations: 147.162 MiB, 38.24% gc time)
```

The julia version runs 7.5 times faster than the matlab version.

## How to Use this

1. [Download latest julia](https://julialang.org/downloads/)
2. start julia. you see something like this:
    ```
    ➜  julia
                   _
       _       _ _(_)_     |  Documentation: https://docs.julialang.org
      (_)     | (_) (_)    |
       _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
      | | | | | | |/ _` |  |
      | | |_| | | | (_| |  |  Version 1.1.0 (2019-01-21)
     _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
    |__/                   |

    julia>
    ```
3. Hit the `]` key to switch to package manager mode. the prompt switches to
    ```
    (v1.1) pkg>
    ```
4. Download this package by pasting this into the `(v1.1) pkg>` prompt and hitting enter.
    ```julia
    dev https://github.com/floswald/DCEGM.jl.git
    ```
5. After this is done, hit backspace or `ctrl-c` to go back to standard `julia>` prompt.
    ```julia
    julia> cd(joinpath(DEPOT_PATH[1],"dev","DCEGM"))  # go to the location of DCEGM on your computer
    ```
6. Go back to package mode by typing `]`. then:
    ```julia
    (v1.1) pkg> activate .     # tell pkg manager to modify current directory
    (DCEGM) pkg> instantiate    # download all dependencies
    (DCEGM) pkg> precompile     # precompile package
    ```
7. Done! :tada: Now try it out. Go back to command mode with `ctrl-c`
    ```julia
    julia> using DCEGM

    julia> m,p = DCEGM.runf();  # runs @fediskhakov version of the algorithm
    ```




## Testing

The package is thoroughly unit tested. Please run `] test` while in the activated project. The main test concerns the file `test/F_test.jl`, where we test the output of this version against the one obtained from @fediskhakov s matlab version up to numerical accuracy. That is, first we save the value and policy functions from the matlab code to ASCII format on disk, then we compute the julia model, then we compare each computed value and policy function.

## Demos

A quick demonstration of how the `upper_env` method works. Given an array of `MLine`s (my version of a *line*, i.e. an array of `x-y` pairs representing a `Point`), this constructs the upper envelope over the lines. Particular attention must be paid to *intersections* between lines.

![](images/demo.png)

Also, points where 2 lines intersect on the initial grid of both lines are *not* intersections.

![](images/demo2.png)
