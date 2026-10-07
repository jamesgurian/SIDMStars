using HDF5

function init_hdf5(filename, n_grid, params)
    h5open(filename, "w") do fid
        # Define initial dimensions and 'unlimited' max dimensions (-1)
        # (n_grid, 0) means we start with zero columns
        initial_dims = (n_grid, 0)
        max_dims = (n_grid, -1)

        # 'chunk' is mandatory for extensible datasets
        # We chunk by 1 column at a time for time-series snapshots
        create_dataset(fid, "t", Float64, ((0,), (-1,)), chunk=(1,))

        for name in ["P", "r", "m", "rho", "s", "L"]
            create_dataset(fid, name, Float64, (initial_dims, max_dims), chunk=(n_grid, 1))
        end

        meta = create_group(fid, "metadata")
        git_hash = try
            readchomp(`git rev-parse HEAD`)
        catch
            "unknown"
        end
        attributes(meta)["git_commit"] = git_hash
        for field in fieldnames(typeof(params))
            val = getfield(params, field)
            attributes(meta)[string(field)] = val isa Union{Number, String} ? val : string(typeof(val))
        end
    end
end

function write_snapshot!(filename, i, t,state)
    h5open(filename, "r+") do fid
        HDF5.set_extent_dims(fid["t"], (i,))
        for name in ["P", "r", "m", "rho", "s", "L"]
            HDF5.set_extent_dims(fid[name], (size(fid[name], 1), i))
        end

        fid["t"][i] = t 
        fid["P"][:, i] = state.P_center  
        fid["r"][:, i] = state.r_center
        fid["m"][:, i] = state.m_center
        fid["rho"][:, i] = state.ρ_center 
        fid["s"][:, i] = state.s_center
        fid["L"][:, i] = 0.5(state.L_face[2:end] .+ state.L_face[1:end-1])
    end
end

"""
    read_snapshot_file(filename)

Read an output file written by [`evolve`](@ref) and return a `Dict`:
- `"t"`: vector of snapshot times.
- `"r"`, `"m"`, `"P"`, `"rho"`, `"s"`, `"L"`: matrices of size
  `(gridpoints - 1, number of snapshots)` at cell centers: radius, enclosed mass,
  pressure, density, entropy `log(P^(3/2)/ρ^(5/2))` and luminosity.
- `"metadata"`: `Dict` of the run's [`RunParams`](@ref) fields.
"""
function read_snapshot_file(filename)
    data = Dict{String, Any}()

    h5open(filename, "r") do fid
        for key in keys(fid)
            obj = fid[key]
            if obj isa HDF5.Group
                data[key] = Dict(k => read_attribute(obj, k) for k in keys(attributes(obj)))
            else
                data[key] = read(obj)
            end
        end
    end

    return data
end