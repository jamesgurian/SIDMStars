function init_hdf5(filename, n_steps, n_grid)
    
    h5open(filename, "w") do fid
    
    # Define the dimensions (Time steps x Spatial grid)
    dims = (n_grid, n_steps)
    
    # Create datasets with Float64 precision
    create_dataset(fid, "t", Float64, (n_steps,), chunk=(n_steps,))
    # We use 'chunking' to optimize how data is written to disk
    create_dataset(fid, "P", Float64, dims, chunk=(n_grid,1))
    create_dataset(fid, "r", Float64, dims, chunk=(n_grid,1))
    create_dataset(fid, "m", Float64, dims, chunk=(n_grid,1))
    create_dataset(fid, "rho", Float64, dims, chunk=(n_grid,1))
    create_dataset(fid, "s", Float64, dims, chunk=(n_grid,1))
    create_dataset(fid, "dsdt", Float64, dims, chunk=(n_grid,1))
    
    end
end

function write_snapshot!(filename, i, t, p, r, m, ρ, s, dsdt)
    # HDF5 allows you to write to specific hyperslabs (slices) on disk
    h5open(filename, "r+") do h5_file
    h5_file["t"][i] = t 
    h5_file["P"][:, i] = p  
    h5_file["r"][:, i] = r
    h5_file["m"][:, i] = m

    h5_file["rho"][:, i] = ρ 
    h5_file["s"][:, i] = s
    h5_file["dsdt"][:, i] = dsdt
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