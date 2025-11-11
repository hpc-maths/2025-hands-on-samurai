import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(14, 14))

# Cell size
dx_coarse = 4.0
dx_fine = dx_coarse / 2.0

# Vertical positions
y_l_plus_1 = 5.0  # Level l+1 (fine) on top
y_l = 2.0  # Level l (coarse) below

cell_height = 1.2

default_fontsize = 18

# ========== Level l+1 (FINE LEVEL - TOP) ==========
ax.text(
    -0.5,
    y_l_plus_1 + cell_height + 0.3,
    "Level l+1 (fine):",
    ha="right",
    fontsize=13,
    weight="bold",
    color="#8B4513",
)

# Real fine cells
fine_cells = [
    (0, y_l_plus_1, dx_fine, cell_height, r"$u_{j-1}^{l+1}$", "#FFD93D", "-", "black"),
    (
        dx_fine,
        y_l_plus_1,
        dx_fine,
        cell_height,
        r"$u_j^{l+1}$",
        "#FFD93D",
        "-",
        "black",
    ),
]

for x, y, w, h, label, color, linestyle, edgecolor in fine_cells:
    rect = patches.Rectangle(
        (x, y),
        w,
        h,
        linewidth=2,
        edgecolor=edgecolor,
        facecolor=color,
        alpha=0.7,
        linestyle=linestyle,
    )
    ax.add_patch(rect)
    ax.text(
        x + w / 2,
        y + h / 2,
        label,
        ha="center",
        va="center",
        fontsize=default_fontsize,
        weight="bold",
    )

# Ghost cell at level l+1 (where level l cell exists)
ghost_x = -dx_fine
rect = patches.Rectangle(
    (ghost_x, y_l_plus_1),
    dx_fine,
    cell_height,
    linewidth=3,
    edgecolor="red",
    facecolor="#FFB6C6",
    alpha=0.5,
    linestyle="--",
)
ax.add_patch(rect)
ax.text(
    ghost_x + dx_fine / 2,
    y_l_plus_1 + cell_height / 2,
    r"$\tilde{u}_{j-2}^{l+1}$",
    ha="center",
    va="center",
    fontsize=default_fontsize,
    weight="bold",
    style="italic",
)
ax.text(
    ghost_x + dx_fine / 2,
    y_l_plus_1 - 0.25,
    "GHOST",
    ha="center",
    va="top",
    fontsize=9,
    weight="bold",
    color="red",
)

# Flux at level l+1 (brown/marron)
flux_x = 0
ax.plot(
    [flux_x, flux_x],
    [y_l_plus_1, y_l_plus_1 + cell_height],
    color="#8B4513",
    linewidth=5,
    linestyle="-",
    zorder=10,
)
ax.text(
    flux_x + 0.15,
    y_l_plus_1 + cell_height / 4,
    r"$F_{j-1/2}^{l+1}$",
    ha="left",
    fontsize=13,
    weight="bold",
    color="#8B4513",
    bbox=dict(
        boxstyle="round,pad=0.3", facecolor="white", edgecolor="#8B4513", linewidth=2
    ),
)

# Arrow showing flux computation from ghost
ax.annotate(
    "",
    xy=(flux_x - 0.1, y_l_plus_1 + cell_height * 0.3),
    xytext=(ghost_x + dx_fine / 2, y_l_plus_1 + cell_height * 0.3),
    arrowprops=dict(arrowstyle="->", color="#8B4513", lw=2.5),
)
ax.text(
    ghost_x + dx_fine * 0.75,
    y_l_plus_1 + cell_height * 0.15,
    r"$F(\tilde{u}_{j-2}^{l+1})$",
    ha="center",
    fontsize=default_fontsize,
    color="#8B4513",
    weight="bold",
)

# Dimension annotation for level l+1
ax.annotate(
    "",
    xy=(dx_fine, y_l_plus_1 - 0.5),
    xytext=(0, y_l_plus_1 - 0.5),
    arrowprops=dict(arrowstyle="<->", color="black", lw=1.5),
)
ax.text(
    dx_fine / 2,
    y_l_plus_1 - 0.7,
    r"$\Delta x_{l+1}$",
    ha="center",
    fontsize=default_fontsize,
    weight="bold",
)

# ========== Level l (COARSE LEVEL - BOTTOM) ==========
ax.text(
    -0.5,
    y_l + cell_height + 0.3,
    "Level l (coarse):",
    ha="right",
    fontsize=13,
    weight="bold",
    color="#2E5090",
)

# Real coarse cell
rect = patches.Rectangle(
    (-dx_coarse, y_l),
    dx_coarse,
    cell_height,
    linewidth=2,
    edgecolor="black",
    facecolor="#A8D5E2",
    alpha=0.7,
    linestyle="-",
)
ax.add_patch(rect)
ax.text(
    -dx_coarse / 2,
    y_l + cell_height / 2,
    r"$u_{i-1}^l$",
    ha="center",
    va="center",
    fontsize=default_fontsize,
    weight="bold",
)

# Ghost cell at level l (where fine cells exist)
rect = patches.Rectangle(
    (0, y_l),
    dx_coarse,
    cell_height,
    linewidth=3,
    edgecolor="red",
    facecolor="#FFB6C6",
    alpha=0.5,
    linestyle="--",
)
ax.add_patch(rect)
ax.text(
    dx_coarse / 2,
    y_l + cell_height / 2,
    r"$\tilde{u}_i^l$",
    ha="center",
    va="center",
    fontsize=default_fontsize,
    weight="bold",
    style="italic",
)
ax.text(
    dx_coarse / 2,
    y_l - 0.25,
    "GHOST",
    ha="center",
    va="top",
    fontsize=9,
    weight="bold",
    color="red",
)

# Flux at level l (blue)
ax.plot(
    [flux_x, flux_x],
    [y_l, y_l + cell_height],
    color="#2E5090",
    linewidth=5,
    linestyle="-",
    zorder=10,
)
ax.text(
    flux_x + 0.15,
    y_l + 3 * cell_height / 4,
    r"$F_{i-1/2}^l$",
    ha="left",
    fontsize=default_fontsize,
    weight="bold",
    color="#2E5090",
    bbox=dict(
        boxstyle="round,pad=0.3", facecolor="white", edgecolor="#2E5090", linewidth=2
    ),
)

# Arrow showing flux computation
ax.annotate(
    "",
    xy=(flux_x - 0.1, y_l + cell_height * 0.7),
    xytext=(-dx_coarse / 2, y_l + cell_height * 0.7),
    arrowprops=dict(arrowstyle="->", color="#2E5090", lw=2.5),
)
ax.text(
    -dx_coarse / 4,
    y_l + cell_height * 0.85,
    r"$F(u_{i-1}^l)$",
    ha="center",
    fontsize=default_fontsize,
    color="#2E5090",
    weight="bold",
)

# Dimension annotation for level l
ax.annotate(
    "",
    xy=(0, y_l - 0.5),
    xytext=(-dx_coarse, y_l - 0.5),
    arrowprops=dict(arrowstyle="<->", color="black", lw=1.5),
)
ax.text(
    -dx_coarse / 2,
    y_l - 0.7,
    r"$\Delta x_l = 2 \Delta x_{l+1}$",
    ha="center",
    fontsize=11,
    weight="bold",
)

# ========== Vertical alignment indicator ==========
# Draw vertical dashed line to show same physical position
ax.plot(
    [flux_x, flux_x],
    [y_l - 0.9, y_l_plus_1 + cell_height + 0.5],
    "r--",
    linewidth=2,
    alpha=0.5,
    zorder=1,
)
# ax.text(
#     flux_x - 0.4,
#     (y_l + y_l_plus_1 + cell_height) / 2,
#     "Same physical\ninterface",
#     ha="center",
#     fontsize=11,
#     color="red",
#     weight="bold",
#     style="italic",
#     bbox=dict(
#         boxstyle="round,pad=0.4", facecolor="white", edgecolor="red", linewidth=2
#     ),
# )

# ========== Explanation box ==========
# problem_text = (
#     r"$\Delta x_l \cdot F(u_{i-1}^l) \neq \Delta x_{l+1} \cdot F(\tilde{u}_{j-2}^{l+1})$"
#     + "\n\n"
#     r"$2\Delta x_{l+1} \cdot F(u_{i-1}^l) \neq \Delta x_{l+1} \cdot F(\tilde{u}_{j-2}^{l+1})$"
#     + "\n\n"
#     "Different ghost values ⇒ Conservation broken!"
# )

# ax.text(
#     3.5,
#     (y_l + y_l_plus_1 + cell_height) / 2,
#     problem_text,
#     ha="left",
#     va="center",
#     fontsize=12,
#     bbox=dict(
#         boxstyle="round,pad=0.8", facecolor="#FFE5E5", edgecolor="red", linewidth=3
#     ),
#     weight="bold",
#     color="darkred",
# )

# # ========== Show what's really there ==========
# # At level l, show the real fine cells
# ax.text(
#     3 * dx_coarse / 4,
#     y_l + cell_height / 4,
#     "Real fine cells\nhere (level l+1)",
#     ha="center",
#     va="center",
#     fontsize=9,
#     style="italic",
#     color="gray",
#     bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
# )

# # At level l+1, show the real coarse cell
# ax.text(
#     ghost_x + dx_fine / 2,
#     y_l + cell_height / 4,
#     "Real coarse\ncell here\n(level l)",
#     ha="center",
#     va="center",
#     fontsize=9,
#     style="italic",
#     color="gray",
#     bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
# )

# Set axis properties
ax.set_xlim(-4, 4)
ax.set_ylim(0.5, y_l_plus_1 + cell_height + 1.2)
ax.set_aspect("equal")
ax.axis("off")

# Add main title
# ax.text(
#     1.5,
#     y_l_plus_1 + cell_height + 0.9,
#     "Why Naive Multi-Resolution Fluxes Break Conservation",
#     ha="center",
#     fontsize=16,
#     weight="bold",
# )

plt.tight_layout()
plt.savefig(
    "flux_conservation_problem.svg",
    bbox_inches="tight",
    dpi=150,
)
print("Diagram saved as: flux_conservation_problem.svg")
