# web_plot.jl
using Bonito
using GLMakie

# Function to create a plot
function open_plot()
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], title="Sine Wave")
    lines!(ax, 0:0.1:10, sin.(0:0.1:10))
    display(fig)
end

# Global reference to track when plot was requested
plot_requested = Ref(false)

# Create the simplest possible app
app = App() do session
    # Simple button that calls a JavaScript function
    DOM.div(
        DOM.h1("Plot Controller"),
        DOM.button(
            "Open Plot",
            id="plotButton",
            style="padding: 10px; background: blue; color: white; border: none; cursor: pointer;"
        ),
        # Basic script to set up the button click handler
        DOM.script("""
            document.getElementById('plotButton').onclick = function() {
                // Add a hidden element to trigger Julia
                const marker = document.createElement('div');
                marker.id = 'plot-requested';
                marker.style.display = 'none';
                document.body.appendChild(marker);
            }
        """)
    )
end

# Start the server
server = Server(app, "127.0.0.1", 8080)

# Main thread checking for the plot request
@async while true
    sleep(0.5)  # Check every half second
    
    # Check if plot was requested by querying for the element
    script = """
        document.getElementById('plot-requested') ? 'true' : 'false'
    """
    result = try
        Bonito.eval_js(server, script)
    catch
        "false"
    end
    
    # If requested and not already shown, show the plot
    if result == "true" && !plot_requested[]
        plot_requested[] = true
    println("Opening plot...")
    open_plot()
        
        # Reset the flag
        Bonito.eval_js(server, """
            const el = document.getElementById('plot-requested');
            if (el) el.parentNode.removeChild(el);
        """)
        plot_requested[] = false
    end
end

println("Server running at http://127.0.0.1:8080")