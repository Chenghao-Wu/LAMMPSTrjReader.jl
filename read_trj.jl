include("/Users/bruce/Dropbox/code/research/hPF-MD/autodiff-hPF-MD/LAMMPSTrjReader/LAMMPSTrjReader.jl")
using .LAMMPSTrjReader
universe        =   LAMMPSTrjReader.Universe("test.lammpstrj")
trj=LAMMPSTrjReader.trajectory(universe,frame_index=1)
println(box(universe))
println(fetch(trj,"atomids"))
println(fetch(trj,"positions"))
println(number_frames(universe))