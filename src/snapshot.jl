using HDF5

function init_hdf5(filename, n_grid)
    h5open(filename, "w") do fid
        # Define initial dimensions and 'unlimited' max dimensions (-1)
        # (n_grid, 0) means we start with zero columns
        initial_dims = (n_grid, 0)
        max_dims = (n_grid, -1) 
        
        # 'chunk' is mandatory for extensible datasets
        # We chunk by 1 column at a time for time-series snapshots
        create_dataset(fid, "t", Float64, ((0,), (-1,)), chunk=(1,))
        
        for name in ["P", "r", "m", "rho", "s"]
            create_dataset(fid, name, Float64, (initial_dims, max_dims), chunk=(n_grid, 1))
        end
    end
end

function write_snapshot!(filename, i, t, p, r, m, ρ, s)
    h5open(filename, "r+") do fid
        HDF5.set_extent_dims(fid["t"], (i,))
        for name in ["P", "r", "m", "rho", "s"]
            HDF5.set_extent_dims(fid[name], (size(fid[name], 1), i))
        end

        fid["t"][i] = t 
        fid["P"][:, i] = p  
        fid["r"][:, i] = r
        fid["m"][:, i] = m
        fid["rho"][:, i] = ρ 
        fid["s"][:, i] = s
    end
end

"""
    read_snapshot_file(filename)
Reads an HDF5 simulation file and returns a Dict of arrays.
"""
function read_snapshot_file(filename)
    data = Dict{String, Any}()
    
    h5open(filename, "r") do fid
        for key in keys(fid)
            data[key] = read(fid[key])
        end
        
    end
    
    return data
end