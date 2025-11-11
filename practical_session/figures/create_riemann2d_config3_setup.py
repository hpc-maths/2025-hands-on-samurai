#!/usr/bin/env python3
"""
Create a schematic diagram of the 2D Riemann problem - Configuration 3 setup.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# Create figure
fig, ax = plt.subplots(figsize=(10, 10))

# Domain boundaries (from riemann_2d.hpp)
domain_width = 1.0
domain_height = 1.0

# Interface positions (from riemann_2d.hpp)
x0 = 0.5
y0 = 0.5

# Draw domain
ax.add_patch(
    Rectangle(
        (0, 0), domain_width, domain_height, fill=False, edgecolor="black", linewidth=2
    )
)

# Draw the four quadrants with different colors
# Quadrant 1: x >= 0.5, y >= 0.5 (top-right)
quad1 = Rectangle(
    (x0, y0),
    domain_width - x0,
    domain_height - y0,
    facecolor="lightblue",
    edgecolor="black",
    linewidth=1.5,
    alpha=0.6,
)
ax.add_patch(quad1)

# Quadrant 2: x < 0.5, y >= 0.5 (top-left)
quad2 = Rectangle(
    (0, y0),
    x0,
    domain_height - y0,
    facecolor="lightcoral",
    edgecolor="black",
    linewidth=1.5,
    alpha=0.6,
)
ax.add_patch(quad2)

# Quadrant 3: x < 0.5, y < 0.5 (bottom-left)
quad3 = Rectangle(
    (0, 0), x0, y0, facecolor="lightgreen", edgecolor="black", linewidth=1.5, alpha=0.6
)
ax.add_patch(quad3)

# Quadrant 4: x >= 0.5, y < 0.5 (bottom-right)
quad4 = Rectangle(
    (x0, 0),
    domain_width - x0,
    y0,
    facecolor="lightyellow",
    edgecolor="black",
    linewidth=1.5,
    alpha=0.6,
)
ax.add_patch(quad4)

# Draw interface lines
ax.plot([x0, x0], [0, domain_height], "k--", linewidth=2, label="Interface")
ax.plot([0, domain_width], [y0, y0], "k--", linewidth=2)

# Add state labels for each quadrant
# Quadrant 1 (top-right)
ax.text(
    0.75,
    0.75,
    r"$\mathbf{Q}_1$"
    + "\n"
    + r"$\rho = 1.5$"
    + "\n"
    + r"$p = 1.5$"
    + "\n"
    + r"$\mathbf{v} = (0, 0)$",
    fontsize=11,
    ha="center",
    va="center",
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", linewidth=1.5),
)

# Quadrant 2 (top-left)
ax.text(
    0.25,
    0.75,
    r"$\mathbf{Q}_2$"
    + "\n"
    + r"$\rho = 0.5323$"
    + "\n"
    + r"$p = 0.3$"
    + "\n"
    + r"$\mathbf{v} = (1.206, 0)$",
    fontsize=11,
    ha="center",
    va="center",
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", linewidth=1.5),
)

# Quadrant 3 (bottom-left)
ax.text(
    0.25,
    0.25,
    r"$\mathbf{Q}_3$"
    + "\n"
    + r"$\rho = 0.138$"
    + "\n"
    + r"$p = 0.29$"
    + "\n"
    + r"$\mathbf{v} = (1.206, 1.206)$",
    fontsize=11,
    ha="center",
    va="center",
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", linewidth=1.5),
)

# Quadrant 4 (bottom-right)
ax.text(
    0.75,
    0.25,
    r"$\mathbf{Q}_4$"
    + "\n"
    + r"$\rho = 0.5323$"
    + "\n"
    + r"$p = 0.3$"
    + "\n"
    + r"$\mathbf{v} = (0, 1.206)$",
    fontsize=11,
    ha="center",
    va="center",
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", linewidth=1.5),
)

# Add velocity arrows to show flow directions
arrow_props = dict(arrowstyle="->", lw=2.5, color="darkblue")

# Q2: velocity (1.206, 0) - rightward (positioned above the text box)
ax.annotate("", xy=(0.45, 0.75), xytext=(0.35, 0.75), arrowprops=arrow_props)

# Q3: velocity (1.206, 1.206) - diagonal (positioned in lower left corner)
ax.annotate("", xy=(0.45, 0.45), xytext=(0.35, 0.35), arrowprops=arrow_props)

# Q4: velocity (0, 1.206) - upward (positioned to the right of the text box)
ax.annotate("", xy=(0.75, 0.45), xytext=(0.75, 0.35), arrowprops=arrow_props)

# Mark interface position
# ax.text(x0 + 0.02, -0.08, r"$x = 0.5$", fontsize=12, ha="left", fontweight="bold")
# ax.text(
#     -0.08, y0, r"$y = 0.5$", fontsize=12, va="center", rotation=90, fontweight="bold"
# )

# Domain labels
# ax.text(domain_width / 2, -0.12, r"$x$", fontsize=14, ha="center", fontweight="bold")
# ax.text(-0.12, domain_height / 2, r"$y$", fontsize=14, va="center", fontweight="bold")

# Add dimensions
ax.annotate(
    "",
    xy=(domain_width, -0.06),
    xytext=(0, -0.06),
    arrowprops=dict(arrowstyle="<->", lw=1.5, color="black"),
)
ax.text(domain_width / 2, -0.09, f"$L = {domain_width}$", fontsize=11, ha="center")

ax.annotate(
    "",
    xy=(-0.06, domain_height),
    xytext=(-0.06, 0),
    arrowprops=dict(arrowstyle="<->", lw=1.5, color="black"),
)
ax.text(
    -0.09,
    domain_height / 2,
    f"$H = {domain_height}$",
    fontsize=11,
    ha="center",
    rotation=90,
    va="center",
)

# Boundary condition label (Neumann on all sides)
bc_text = "Boundary conditions: Neumann on all sides"
ax.text(
    0.5,
    1.08,
    bc_text,
    fontsize=11,
    ha="center",
    style="italic",
    bbox=dict(boxstyle="round,pad=0.5", facecolor="lavender", alpha=0.8),
)

# Set axis properties
ax.set_xlim(-0.15, domain_width + 0.05)
ax.set_ylim(-0.22, domain_height + 0.12)
ax.set_aspect("equal")
ax.axis("off")

plt.tight_layout()
plt.savefig("riemann2d_config3_setup.svg", bbox_inches="tight")
print("Figure saved: riemann2d_config3_setup.svg")
plt.close()

print("2D Riemann Configuration 3 setup diagram created successfully!")
