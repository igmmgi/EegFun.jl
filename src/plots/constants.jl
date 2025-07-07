# Plotting constants for eegfun
# This file consolidates all default constants used across plotting functions

# Trigger overview plot constants
const DEFAULT_MARKER_SIZE = 15
const DEFAULT_LINE_WIDTH_TRIGGER = 1
const DEFAULT_LINE_OFFSET = 0.1

# Layout plot constants
const DEFAULT_HEAD_RADIUS = 88.0  # mm - standard head size for EEG layouts
const DEFAULT_BASE_SIZE = 15      # pixels - default electrode point size
const DEFAULT_HOVER_SIZE = 25     # pixels - electrode size on hover
const DEFAULT_LINE_WIDTH_LAYOUT = 3  # pixels - width of connection lines
const DEFAULT_HEAD_COLOR = :black
const DEFAULT_HEAD_LINEWIDTH = 2
const DEFAULT_POINT_COLOR = :black
const DEFAULT_POINT_MARKER = :circle
const DEFAULT_POINT_SIZE = 12
const DEFAULT_LABEL_COLOR = :black
const DEFAULT_LABEL_FONTSIZE = 20

# Trigger timing plot constants
const DEFAULT_TIMELINE_WIDTH = 2
const DEFAULT_EVENT_LINE_WIDTH = 1
const DEFAULT_FONT_SIZE = 24
const DEFAULT_WINDOW_SIZE = 10.0  # seconds

# Layout-specific constants
const HEAD_EAR_RATIO = 1/7        # Ratio of ear size to head radius
const HEAD_NOSE_SCALE = 4.0       # Scale factor for nose size
const MAX_REASONABLE_COORDINATE = 1000.0  # mm - maximum reasonable coordinate value

# Trigger timing specific constants
const EVENT_LINE_HEIGHT = 0.05
const TIME_LABEL_OFFSET = -0.001
const Y_MIN_LIMIT = -0.1
const Y_MAX_LIMIT = 0.15
const MIN_WINDOW_SIZE = 1.0       # seconds
const MAX_WINDOW_SIZE = 60.0      # seconds
const WINDOW_SIZE_STEP = 1.0      # seconds
const POSITION_STEP = 0.5         # seconds 