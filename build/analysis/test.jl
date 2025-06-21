using JLD2, CodecZstd, Parquet2, Arrow, DataFrames

# Create data that should compress well (repeated patterns)
test_data = repeat([1.0, 2.0, 3.0, 4.0, 5.0], 1000, 100)  # Lots of repetition

# Save both versions
jldsave("test_uncompressed.jld2"; data=test_data)
jldsave("test_compressed.jld2"; data=test_data, compress=ZstdCompressor(level=10))

# Compare sizes
println("Uncompressed: ", filesize("test_uncompressed.jld2"), " bytes")
println("Compressed: ", filesize("test_compressed.jld2"), " bytes")

# Create some test EEG-like data
function create_test_eeg_data(n_channels=64, n_samples=1000)
    # Simulate EEG data with channel names and some realistic patterns
    channels = ["CH$i" for i in 1:n_channels]
    
    # Create semi-realistic EEG data (not completely random)
    data = Dict()
    for ch in channels
        # EEG-like signal: low frequency base + some noise
        t = 1:n_samples
        signal = 0.1 * sin.(2π * 0.01 * t) + 0.05 * randn(n_samples)
        data[ch] = signal
    end
    
    # Add some metadata columns
    data["time"] = collect(1:n_samples) ./ 1000.0  # Time in seconds
    data["is_artifact"] = rand(Bool, n_samples)
    
    return DataFrame(data)
end

# Test function
function compare_formats()
    # Create test data
    test_df = create_test_eeg_data()
    
    println("Original DataFrame size in memory: ", Base.summarysize(test_df), " bytes")
    
    # Test JLD2
    jldsave("test_data.jld2"; data=test_df)
    jld2_size = filesize("test_data.jld2")
    println("JLD2 file size: ", jld2_size, " bytes")
    
    # Test Parquet
    Parquet2.writefile("test_data.parquet", test_df)
    parquet_size = filesize("test_data.parquet")
    println("Parquet file size: ", parquet_size, " bytes")
    
    # Test Arrow
    Arrow.write("test_data.arrow", test_df)
    arrow_size = filesize("test_data.arrow")
    println("Arrow file size: ", arrow_size, " bytes")
    
    # Calculate compression ratios
    println("\nCompression ratios (smaller is better):")
    println("JLD2: ", round(jld2_size / Base.summarysize(test_df), digits=3))
    println("Parquet: ", round(parquet_size / Base.summarysize(test_df), digits=3))
    println("Arrow: ", round(arrow_size / Base.summarysize(test_df), digits=3))
    
    # Clean up
    rm("test_data.jld2")
    rm("test_data.parquet")
    rm("test_data.arrow")
end

# Run the test
compare_formats()