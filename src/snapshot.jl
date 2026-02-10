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
        
        # Grid-based datasets
        for name in ["P", "r", "m", "rho", "s", "dsdt"]
            create_dataset(fid, name, Float64, (initial_dims, max_dims), chunk=(n_grid, 1))
        end
    end
end

function write_snapshot!(filename, i, t, p, r, m, ρ, s, dsdt)
    h5open(filename, "r+") do fid
        # 1. Resize datasets to accommodate the i-th entry
        HDF5.set_extent_dims(fid["t"], (i,))
        for name in ["P", "r", "m", "rho", "s", "dsdt"]
            HDF5.set_extent_dims(fid[name], (size(fid[name], 1), i))
        end

        # 2. Write the data to the newly created column
        fid["t"][i] = t 
        fid["P"][:, i] = p  
        fid["r"][:, i] = r
        fid["m"][:, i] = m
        fid["rho"][:, i] = ρ 
        fid["s"][:, i] = s
        fid["dsdt"][:, i] = dsdt
    end
end

"""
    read_snapshot_file(filename)
Reads an HDF5 simulation file and returns a Dict of arrays.
"""
function read_snapshot_file(filename)
    data = Dict{String, Any}()
    
    h5open(filename, "r") do fid
        # Iterate through every dataset in the file
        for key in keys(fid)
            # read(fid[key]) loads the entire dataset into a Julia Array
            data[key] = read(fid[key])
        end
        
        # Optional: Read metadata attributes if you saved them
        # data["description"] = read_attribute(fid, "description")
    end
    
    return data
end