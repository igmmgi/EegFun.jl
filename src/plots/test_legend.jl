using GLMakie
using eegfun

# Create a figure and axis
fig = Figure()
ax = Axis(fig[1, 1])

# Create some test plots with labels
lines!(ax, 1:10, 1:10, label = "Line 1")
lines!(ax, 1:10, 2:11, label = "Line 2")
lines!(ax, 1:10, 3:12, label = "Line 3")

# Test 1: Basic legend with position
println("Test 1: Basic legend with position")
axislegend(ax; position = :lt)
display(fig)

# Test 2: Legend with position and framevisible
println("\nTest 2: Legend with position and framevisible")
fig2 = Figure()
ax2 = Axis(fig2[1, 1])
lines!(ax2, 1:10, 1:10, label = "Line 1")
lines!(ax2, 1:10, 2:11, label = "Line 2")
# Simulate what plot_erp does
plot_kwargs = Dict(
    :legend => true,
    :legend_position => (0.5, 0.5),
    :legend_framevisible => true,
    :legend_label => "Test Legend"
)
legend_kwargs = eegfun._extract_legend_kwargs!(plot_kwargs)
println("Extracted legend_kwargs: ", legend_kwargs)
println("Remaining plot_kwargs: ", plot_kwargs)
axislegend(ax2, "help"; position = plot_kwargs[:legend_position], legend_kwargs...)
display(fig2)

# Test 3: Test with multiple legend attributes
println("\nTest 3: Multiple legend attributes")
fig3 = Figure()
ax3 = Axis(fig3[1, 1])
lines!(ax3, 1:10, 1:10, label = "Line 1")
lines!(ax3, 1:10, 2:11, label = "Line 2")

plot_kwargs3 = Dict(
    :legend => true,
    :legend_position => (:left, :bottom),
    :legend_framevisible => true,
    :legend_bgcolor => :lightgray,
    :legend_padding => (10, 10, 10, 10),
    :legend_label => "Custom Legend"
)

legend_kwargs3 = eegfun._extract_legend_kwargs!(plot_kwargs3)
println("Extracted legend_kwargs: ", legend_kwargs3)

axislegend(ax3, plot_kwargs3[:legend_label]; position = plot_kwargs3[:legend_position], legend_kwargs3...)
display(fig3)

println("\nAll tests completed!")
