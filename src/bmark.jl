using DCEGM
using BenchmarkTools
println("using $(Threads.nthreads()) threads")
println("===========================")
na = 5000; ny = 50

println("fedors model on na=$na,ny=$ny")
DCEGM.runf(par = Dict(:na => na, :ny => ny));
@btime DCEGM.runf(par = Dict(:na => na, :ny => ny));

na = 500; ny = 15
println("general model on na=$na,ny=$ny")
DCEGM.rung(par = Dict(:na => na, :ny => ny));
@btime DCEGM.rung(par = Dict(:na => na, :ny => ny));
