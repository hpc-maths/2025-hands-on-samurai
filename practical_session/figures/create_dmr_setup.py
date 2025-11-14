#!/usr/bin/env python3
"""
Create a schematic diagram of the Double Mach Reflection problem setup.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch
from matplotlib.patches import Arc
import matplotlib.patches as mpatches

# Create figure
fig, ax = plt.subplots(figsize=(14, 7))

# Domain boundaries
domain_width = 4.0
domain_height = 1.0

# Draw domain
# ax.add_patch(
#     Rectangle(
#         (0, 0), domain_width, domain_height, fill=False, edgecolor="black", linewidth=2
#     )
# )

# Shock angle: 60 degrees from horizontal (pi/3 from code)
shock_angle = np.pi / 3  # 60 degrees in radians

# Initial shock position at x = 2/3 (from code)
x_shock_init = 2.0 / 3.0

# Shock line extends from point (1/6, 0) at 30 degrees
shock_length = 2.0
x_shock_end = x_shock_init + shock_length * np.cos(shock_angle)
y_shock_end = shock_length * np.sin(shock_angle)

# Draw the shock wave
shock_line = ax.plot(
    [x_shock_init, x_shock_end],
    [0, y_shock_end],
    "r-",
    linewidth=3,
    label="Initial shock position",
)

# Draw BC
x1 = x_shock_init + 1 / np.tan(shock_angle)
ax.plot(
    [0, x1],
    [domain_height, domain_height],
    "r-",
    linewidth=2,
)  # Top
ax.plot(
    [x1, domain_width],
    [domain_height, domain_height],
    "b-",
    linewidth=2,
)  # Top
ax.plot([domain_width, domain_width], [0, domain_height], "g-", linewidth=2)  # Right
ax.plot([0, 0], [0, domain_height], "r-", linewidth=2)  # Left
ax.plot([0, x_shock_init], [0, 0], "r-", linewidth=2)  # Bottom
ax.plot([x_shock_init, domain_width], [0, 0], "k-", linewidth=2)  # Bottom

# Add arrow to show shock direction
arrow_start_x = x_shock_init + 0.5 * np.cos(shock_angle)
arrow_start_y = 0.5 * np.sin(shock_angle)
arrow_dx = 0.3 * np.cos(shock_angle)
arrow_dy = 0.3 * np.sin(shock_angle)

ax.arrow(
    arrow_start_x,
    arrow_start_y,
    arrow_dx,
    arrow_dy,
    head_width=0.08,
    head_length=0.1,
    fc="red",
    ec="red",
    linewidth=2,
)

# Draw reflecting wall (from x=0 to x=2/3)
wall_start = 2.0 / 3.0
wall_end = 4.0
ax.plot([wall_start, wall_end], [0, 0], "k-", linewidth=5)

# Add hatching below the wall to show it's solid
hatch_depth = -0.1
ax.fill_between(
    [wall_start, wall_end],
    [0, 0],
    [hatch_depth, hatch_depth],
    hatch="///",
    facecolor="lightgray",
    edgecolor="black",
)

# Add angle arc
angle_radius = 0.4
angle_arc = Arc(
    (x_shock_init, 0),
    2 * angle_radius,
    2 * angle_radius,
    angle=0,
    theta1=0,
    theta2=60,
    color="blue",
    linewidth=2,
)
ax.add_patch(angle_arc)
ax.text(
    x_shock_init + 0.15, 0.12, r"$60°$", fontsize=14, color="blue", fontweight="bold"
)

# Label regions
# Pre-shock state (right of shock)
ax.text(
    2.8,
    0.75,
    "Pre-shock state\n" + r"$(\rho, u, v, p) = (1.4, 0, 0, 1)$",
    fontsize=11,
    bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8),
    ha="center",
)

# Post-shock state (left of shock, above wall)
# From code: rho=8, p=116.5, v=(8.25*sin(60°), -8.25*cos(60°))
post_shock_x = 0.5
post_shock_y = 0.75
ax.text(
    post_shock_x,
    post_shock_y,
    "Post-shock state\n"
    + r"$(\rho, p) = (8, 116.5)$"
    + "\n"
    + r"$\mathbf{v} = (8.25\sin 60°, -8.25\cos 60°)$",
    fontsize=11,
    bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.8),
    ha="center",
)

# Add Mach number label on shock
ax.text(
    x_shock_init + 0.7,
    0.35,
    "M = 10",
    fontsize=13,
    color="red",
    fontweight="bold",
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="red", linewidth=2),
)

# Domain labels
# ax.text(domain_width / 2, -0.25, r"$x$", fontsize=14, ha="center", fontweight="bold")
# ax.text(-0.15, domain_height / 2, r"$y$", fontsize=14, va="center", fontweight="bold")

# Add dimensions
ax.annotate(
    "",
    xy=(domain_width, -0.2),
    xytext=(0, -0.2),
    arrowprops=dict(arrowstyle="<->", lw=1.5, color="black"),
)
ax.text(domain_width / 2, -0.27, f"L = {domain_width}", fontsize=11, ha="center")

ax.annotate(
    "",
    xy=(-0.1, domain_height),
    xytext=(-0.1, 0),
    arrowprops=dict(arrowstyle="<->", lw=1.5, color="black"),
)
ax.text(
    -0.2,
    domain_height / 2,
    f"H = {domain_height}",
    fontsize=11,
    ha="center",
    rotation=90,
    va="center",
)

# Mark initial shock position
ax.plot([x_shock_init, x_shock_init], [0, -0.05], "k--", linewidth=1)
ax.text(x_shock_init, -0.15, r"$x = 2/3$", fontsize=10, ha="center")

# Boundary condition labels
# Top boundary
# ax.text(
#     domain_width / 2,
#     domain_height + 0.12,
#     "Inflow (post-shock)",
#     fontsize=10,
#     ha="center",
#     style="italic",
#     bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.6),
# )

# Right boundary
ax.text(
    domain_width + 0.15,
    domain_height / 2,
    "Outflow",
    fontsize=10,
    va="center",
    rotation=270,
    style="italic",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightcoral", alpha=0.6),
)

# Bottom boundary (right of wall)
ax.text(
    (wall_start + domain_width) / 2,
    0.1,
    "Reflecting wall",
    fontsize=10,
    ha="center",
    style="italic",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="gray", alpha=0.6),
)

# Set axis properties
ax.set_xlim(-0.5, domain_width + 0.5)
ax.set_ylim(-0.6, domain_height + 0.25)
ax.set_aspect("equal")
ax.axis("off")

# Add legend
# ax.legend(loc="upper right", fontsize=11, framealpha=0.9)

# Add note
note_text = (
    "Note: Shock moves from left to right.\n"
    "Complex reflection pattern forms at the wall."
)
ax.text(
    0.02,
    0.98,
    note_text,
    transform=ax.transAxes,
    fontsize=9,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7),
)

plt.tight_layout()
plt.savefig("double_mach_reflection_setup.svg", bbox_inches="tight")
print("Figure saved: double_mach_reflection_setup.svg")
plt.close()

print("Double Mach Reflection setup diagram created successfully!")
