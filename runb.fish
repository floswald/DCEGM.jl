 #! /usr/local/bin/fish
for n in (seq 2)
    env JULIA_NUM_THREADS=$n julia --project=. ./src/bmark.jl   
end
