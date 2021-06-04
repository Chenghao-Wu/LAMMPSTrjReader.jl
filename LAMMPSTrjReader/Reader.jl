export  Universe,
        box,
        trajectory,
        fetch,
        number_frames

# Sections will all start with one of these words
# and run until the next section title
SECTIONS = [
    "ITEM: TIMESTEP",  # Molecular topology sections
]

HEADERS = [
    "type",
    "id",
    "x",
    "y",
    "z",
    "ix",
    "iy",
    "iz",
    "vx",
    "vy",
    "vz",
]

function get_StartLineNumber(lines)
    starts=[]
    Sects=Dict{Int64,Any}()
    for (line_index,line) in enumerate(lines)
        line=split(line,"#")
        if length(line)>0
            if strip(line[1]) in SECTIONS
                #println(line[1])
                push!(starts,line_index)
            end
        end
    end
    push!(starts,length(lines))
    #print(starts)
    for (i,sectii) in enumerate(starts[1:length(starts)-1])
        push!(Sects,i=>lines[sectii:starts[i+1]-1])
    end
    #print(Sects[1])
    return Sects
end


function _parse_box(sects,frame_index)
    x1, x2 = (parse(Float64,x) for x in split(sects[frame_index][6]))
    x = x2 - x1
    y1, y2 = (parse(Float64,x) for x in split(sects[frame_index][7]))
    y = y2 - y1
    z1, z2 = (parse(Float64,x) for x in split(sects[frame_index][8]))
    z = z2 - z1

    unitcell = zeros(3)
    unitcell[1:3] .= x, y, z
    return unitcell
end

function _parse_trajectory(sects,headers)
    number_atoms=length(sects)-9
    #print(sects[10:number_atoms+9])
    atom_ids = zeros(Int64,number_atoms)
    atom_types = zeros(Int64,number_atoms)
    positions=zeros(3,number_atoms)
    images=zeros(3,number_atoms)
    velocities=zeros(3,number_atoms)

    for (line_index,line) in enumerate(sects[10:number_atoms+9])
        line=split(line)
        for (item_index,item) in enumerate(line)
            if headers[item_index]=="id"
                atom_ids[line_index]=parse(Int64,strip(item))
            end
            if headers[item_index]=="type"
                atom_types[line_index]=parse(Int64,strip(item))
            end
            if headers[item_index]=="x"
                positions[1,line_index]=parse(Float64,strip(item))
            end
            if headers[item_index]=="y"
                positions[2,line_index]=parse(Float64,strip(item))
            end
            if headers[item_index]=="z"
                positions[3,line_index]=parse(Float64,strip(item))
            end
            if headers[item_index]=="ix"
                images[1,line_index]=parse(Float64,strip(item))
            end
            if headers[item_index]=="iy"
                images[2,line_index]=parse(Float64,strip(item))
            end
            if headers[item_index]=="iz"
                images[3,line_index]=parse(Float64,strip(item))
            end
            if headers[item_index]=="vx"
                velocities[1,line_index]=parse(Float64,strip(item))
            end
            if headers[item_index]=="vy"
                velocities[2,line_index]=parse(Float64,strip(item))
            end
            if headers[item_index]=="vz"
                velocities[3,line_index]=parse(Float64,strip(item))
            end
        end
    end

    order = sortperm(atom_ids)
    atomids=atom_ids[order]
    atomtypes=atom_types[order]
    positions_order=positions[order]
    images_order=images[order]
    velocities_order=velocities[order]
    
    return atomids,atomtypes,positions_order,images_order,velocities_order
end

function get_TokenLineTypes(sects)
    headers=Dict{Int64,Any}()
    TokenLineTypes=split(sects[1][9])
    for (line_index,line) in enumerate(TokenLineTypes)
        #line=split(line,"#")
        if length(line)>0
            for token in HEADERS
                if strip(line) == token
                    push!(headers,line_index-2=>line)
                end
            end
        end
    end
    return headers
end

struct Universe
    sects
    headers
end


function Universe(filename::String)
    file_in=open(filename)
    lines=readlines(file_in)
    sects=get_StartLineNumber(lines)
    headers=get_TokenLineTypes(sects)
    #print(headers)
    #sections=keys(sects)
    #string=""
    #for sectii in sections
    #    string*=sectii*" "
    #end
    #@info "Obtaining Sections "*string
    return Universe(sects,headers)
end

#try
#    masses=_parse_mass(sects["Masses"])
#catch e
#    @error "LAMMPSTrjReader.jl can not find Masses for atoms"
#end
#
#if "Atoms" ∉ keys(sects)
#    @error "LAMMPSTrjReader.jl can not find Atoms"
#end
#
#try 
#    atoms=_parse_atoms(sects["Atoms"],false,masses)
#catch e
#    @error "LAMMPSTrjReader.jl Failed to parse atoms section."
#end

function box(universe;frame_index=1)
    box=nothing
    try 
        box=_parse_box(universe.sects,frame_index)
    catch e
        @error "LAMMPSTrjReader.jl Failed to parse box header at Frame "*string(frame_index)
    end
    box
end

trajectory_order=Dict{String,Any}("atomids"=>1,"atomtypes"=>2,"positions"=>3,"images"=>4,"velocities"=>5)

function trajectory(universe;frame_index=1)
    try #atomids,atomtypes,positions_order,images_order,velocities_order
        trajectories=_parse_trajectory(universe.sects[frame_index],universe.headers)
        return trajectories
    catch e
        @error "LAMMPSTrjReader.jl Failed to parse trajectory at Frame "*string(frame_index)
    end
end

function fetch(trajectories,item)
    index=trajectory_order[item]
    return trajectories[index]
end


function number_frames(universe)
    return length(keys(universe.sects))
end