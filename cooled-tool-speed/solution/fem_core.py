"""Steady plane-strain heat conduction and thermal stress for a cooled tool tip."""

from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.spatial import Delaunay


def conductivity(temperature_c, spec):
    rise = np.maximum(np.asarray(temperature_c, dtype=float) - spec["conductivity_ref_temperature_c"], 0.0)
    return spec["conductivity_ref_w_mk"] / (1.0 + spec["conductivity_softening_per_k"] * rise)


def _level_counts(level):
    if level == "medium":
        return {"n_side": 8, "n_rad": 5, "n_left": 8, "n_right": 10, "n_bottom": 4, "n_top": 6}
    if level == "fine":
        return {"n_side": 16, "n_rad": 10, "n_left": 16, "n_right": 20, "n_bottom": 8, "n_top": 12}
    raise ValueError(level)


def build_mesh(spec, level):
    if spec["hole_radius_m"] <= 0.0:
        return _cartesian_mesh(spec, level)
    counts = _level_counts(level)
    half = 0.0023
    center = np.array([spec["hole_center_x_m"], spec["hole_center_y_m"]])
    square = np.array([center[0] - half, center[0] + half, center[1] - half, center[1] + half])
    welder = _NodeWelder()
    circle_ids = []
    quads = []
    n_side = counts["n_side"]
    n_rad = counts["n_rad"]
    perimeter = _square_perimeter(square, n_side)
    ring = []
    for point in perimeter:
        direction = point - center
        direction = direction / np.linalg.norm(direction)
        column = []
        for step in range(n_rad + 1):
            blend = step / n_rad
            xy = center + ((1.0 - blend) * spec["hole_radius_m"] + blend * np.linalg.norm(point - center)) * direction
            column.append(welder.add(*xy))
        circle_ids.append(column[0])
        ring.append(column)
    for index in range(len(ring)):
        nxt = ring[(index + 1) % len(ring)]
        for step in range(n_rad):
            quads.append([ring[index][step], nxt[step], nxt[step + 1], ring[index][step + 1]])
    blocks = _outer_blocks(spec, square, counts)
    for block in blocks:
        _add_cartesian_block(welder, quads, *block)
    points = welder.points()
    triangles = _quads_to_triangles(points, quads)
    return _boundaries(points, triangles, np.array(circle_ids, dtype=int), spec)


def _cartesian_mesh(spec, level):
    counts = _level_counts(level)
    nx = counts["n_left"] + counts["n_side"] + counts["n_right"]
    ny = counts["n_bottom"] + counts["n_side"] + counts["n_top"]
    if level == "fine":
        nx *= 1
        ny *= 1
    welder = _NodeWelder()
    quads = []
    _add_cartesian_block(
        welder,
        quads,
        0.0,
        spec["width_m"],
        0.0,
        spec["height_m"],
        nx,
        ny,
    )
    triangles = _quads_to_triangles(welder.points(), quads)
    return _boundaries(welder.points(), triangles, np.array([], dtype=int), spec)


class _NodeWelder:
    def __init__(self):
        self._index = {}
        self._points = []

    def add(self, x, y):
        key = (round(float(x), 10), round(float(y), 10))
        if key not in self._index:
            self._index[key] = len(self._points)
            self._points.append([key[0], key[1]])
        return self._index[key]

    def points(self):
        return np.array(self._points, dtype=float)


def _square_perimeter(square, n_side):
    x0, x1, y0, y1 = square
    points = []
    for y in np.linspace(y0, y1, n_side, endpoint=False):
        points.append(np.array([x1, y]))
    for x in np.linspace(x1, x0, n_side, endpoint=False):
        points.append(np.array([x, y1]))
    for y in np.linspace(y1, y0, n_side, endpoint=False):
        points.append(np.array([x0, y]))
    for x in np.linspace(x0, x1, n_side, endpoint=False):
        points.append(np.array([x, y0]))
    return points


def _outer_blocks(spec, square, counts):
    x0, x1, y0, y1 = square
    width = spec["width_m"]
    height = spec["height_m"]
    return [
        (x1, width, y0, y1, counts["n_right"], counts["n_side"]),
        (0.0, x0, y0, y1, counts["n_left"], counts["n_side"]),
        (x0, x1, y1, height, counts["n_side"], counts["n_top"]),
        (x0, x1, 0.0, y0, counts["n_side"], counts["n_bottom"]),
        (0.0, x0, y1, height, counts["n_left"], counts["n_top"]),
        (x1, width, y1, height, counts["n_right"], counts["n_top"]),
        (0.0, x0, 0.0, y0, counts["n_left"], counts["n_bottom"]),
        (x1, width, 0.0, y0, counts["n_right"], counts["n_bottom"]),
    ]


def _add_cartesian_block(welder, quads, x0, x1, y0, y1, nx, ny):
    ids = [[welder.add(x, y) for x in np.linspace(x0, x1, nx + 1)] for y in np.linspace(y0, y1, ny + 1)]
    for j in range(ny):
        for i in range(nx):
            quads.append([ids[j][i], ids[j][i + 1], ids[j + 1][i + 1], ids[j + 1][i]])


def _quads_to_triangles(points, quads):
    triangles = []
    for a, b, c, d in quads:
        triangles.append([a, b, c])
        triangles.append([a, c, d])
    triangles = np.array(triangles, dtype=int)
    signed = _signed_area(points, triangles)
    triangles[signed < 0.0] = triangles[signed < 0.0][:, ::-1]
    return triangles[np.abs(signed) > 1e-16]


def _signed_area(points, triangles):
    a = points[triangles[:, 0]]
    b = points[triangles[:, 1]]
    c = points[triangles[:, 2]]
    return 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))


def _boundaries(points, triangles, circle_ids, spec):
    edge_count = {}
    for tri in triangles:
        for i, j in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            key = (int(i), int(j)) if i < j else (int(j), int(i))
            edge_count[key] = edge_count.get(key, 0) + 1
    boundary = [key for key, count in edge_count.items() if count == 1]
    circle = set(int(node) for node in circle_ids)
    width = spec["width_m"]
    height = spec["height_m"]
    patch = spec["heat_patch_length_m"]
    hole_edges = []
    flux_edges = []
    air_edges = []
    for i, j in boundary:
        if i in circle and j in circle:
            hole_edges.append((i, j))
            continue
        midpoint = 0.5 * (points[i] + points[j])
        on_shank = midpoint[0] >= width - 1e-9
        if on_shank:
            continue
        on_rake = midpoint[1] >= height - 1e-9 and midpoint[0] <= patch + 1e-9
        if on_rake:
            flux_edges.append((i, j))
        else:
            air_edges.append((i, j))
    shank_nodes = np.where(points[:, 0] >= width - 1e-9)[0]
    return {
        "points": points,
        "triangles": triangles,
        "hole_edges": np.array(hole_edges, dtype=int).reshape(-1, 2),
        "flux_edges": np.array(flux_edges, dtype=int).reshape(-1, 2),
        "air_edges": np.array(air_edges, dtype=int).reshape(-1, 2),
        "shank_nodes": shank_nodes.astype(int),
    }


def _apply_dirichlet(matrix, rhs, nodes, values, stride):
    dense_nodes = []
    for node, value in zip(nodes, values):
        for offset in range(stride):
            dense_nodes.append(node * stride + offset)
    dense_nodes = np.array(dense_nodes, dtype=int)
    prescribed = np.zeros(matrix.shape[0])
    cursor = 0
    for value in values:
        for offset in range(stride):
            prescribed[dense_nodes[cursor]] = value if stride == 1 else value[offset]
            cursor += 1
    rhs -= matrix @ prescribed
    rhs[dense_nodes] = prescribed[dense_nodes]
    matrix = matrix.tolil()
    matrix[dense_nodes, :] = 0.0
    matrix[:, dense_nodes] = 0.0
    matrix[dense_nodes, dense_nodes] = 1.0
    return matrix.tocsr(), rhs


def solve_temperature(mesh, spec, speed_m_s):
    points = mesh["points"]
    triangles = mesh["triangles"]
    count = len(points)
    temperature = np.full(count, spec["stress_free_temperature_c"])
    flux = spec["heat_flux_per_speed_w_m2_per_m_s"] * speed_m_s
    for _ in range(16):
        rows, cols, vals = [], [], []
        rhs = np.zeros(count)
        element_k = conductivity(temperature[triangles].mean(axis=1), spec)
        for tri, conduct in zip(triangles, element_k):
            local, area = _conduction(points[tri], conduct)
            _scatter(rows, cols, vals, tri, local, 1)
        _robin(points, mesh["hole_edges"], spec["coolant_h_w_m2k"], spec["coolant_temperature_c"], rows, cols, vals, rhs)
        _robin(points, mesh["air_edges"], spec["air_h_w_m2k"], spec["air_temperature_c"], rows, cols, vals, rhs)
        _neumann(points, mesh["flux_edges"], flux, rhs)
        matrix = coo_matrix((vals, (rows, cols)), shape=(count, count)).tocsr()
        shank = mesh["shank_nodes"]
        matrix, rhs = _apply_dirichlet(matrix, rhs, shank, np.full(len(shank), spec["shank_temperature_c"]), 1)
        updated = spsolve(matrix, rhs)
        updated = 0.6 * updated + 0.4 * temperature
        if np.max(np.abs(updated - temperature)) < 0.02:
            temperature = updated
            break
        temperature = updated
    return temperature


def solve_displacement(mesh, spec, temperature):
    points = mesh["points"]
    triangles = mesh["triangles"]
    count = len(points)
    ndof = 2 * count
    modulus = _plane_strain_modulus(spec["young_modulus_pa"], spec["poisson_ratio"])
    thermal_scale = spec["alpha_per_c"] * (1.0 + spec["poisson_ratio"])
    rows, cols, vals = [], [], []
    rhs = np.zeros(ndof)
    for tri in triangles:
        coords = points[tri]
        area, grads = _gradients(coords)
        body = modulus @ np.array([1.0, 1.0, 0.0])
        mean_rise = temperature[tri].mean() - spec["stress_free_temperature_c"]
        local = np.zeros((6, 6))
        force = np.zeros(6)
        b_matrix = _strain_matrix(grads)
        local += area * (b_matrix.T @ modulus @ b_matrix)
        force += area * thermal_scale * mean_rise * (b_matrix.T @ body)
        dofs = np.array([2 * tri[0], 2 * tri[0] + 1, 2 * tri[1], 2 * tri[1] + 1, 2 * tri[2], 2 * tri[2] + 1])
        _scatter(rows, cols, vals, dofs, local, 1)
        rhs[dofs] += force
    matrix = coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr()
    shank = mesh["shank_nodes"]
    zeros = np.zeros((len(shank), 2))
    matrix, rhs = _apply_dirichlet(matrix, rhs, shank, zeros, 2)
    displacement = spsolve(matrix, rhs).reshape(-1, 2)
    return displacement


def sample_rake_temperatures(mesh, temperature, spec):
    fractions = (0.15, 0.35, 0.55, 0.75, 0.95)
    values = []
    height = spec["height_m"]
    for fraction in fractions:
        target = np.array([fraction * spec["heat_patch_length_m"], height])
        values.append(_sample_top(mesh, temperature, target))
    return values


def cutting_corner_growth_um(mesh, displacement):
    points = mesh["points"]
    corner = int(np.argmin(points[:, 0] ** 2 + (points[:, 1] - points[:, 1].max()) ** 2))
    return float(-displacement[corner, 0] * 1e6)


def limit_speed(mesh, spec):
    limit = spec["edge_temperature_limit_c"]
    corner = _corner_index(mesh)

    def corner_temperature(speed):
        temperature = solve_temperature(mesh, spec, speed)
        return float(temperature[corner]), temperature

    low_temperature, _ = corner_temperature(spec["speed_bracket_low_m_s"])
    high_temperature, _ = corner_temperature(spec["speed_bracket_high_m_s"])
    if not (low_temperature < limit < high_temperature):
        raise RuntimeError(f"limit {limit} is outside bracket temperatures {low_temperature}, {high_temperature}")
    low = spec["speed_bracket_low_m_s"]
    high = spec["speed_bracket_high_m_s"]
    for _ in range(28):
        mid = 0.5 * (low + high)
        value, _ = corner_temperature(mid)
        if value < limit:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def evaluate(spec, level, with_hole=True, constant_conductivity=False):
    task = dict(spec)
    if not with_hole:
        task["hole_radius_m"] = 0.0
    if constant_conductivity:
        task["conductivity_softening_per_k"] = 0.0
    mesh = build_mesh(task, level)
    temperature = solve_temperature(mesh, task, task["reference_speed_m_s"])
    displacement = solve_displacement(mesh, task, temperature)
    return {
        "rake_temperature_c": sample_rake_temperatures(mesh, temperature, task),
        "edge_growth_um": cutting_corner_growth_um(mesh, displacement),
        "limit_speed_m_s": limit_speed(mesh, task),
        "corner_temperature_c": float(temperature[_corner_index(mesh)]),
        "node_count": int(len(mesh["points"])),
        "element_count": int(len(mesh["triangles"])),
    }


def annulus_temperature_error_c():
    """Constant-k annular wall: inner convection, outer prescribed temperature."""
    inner = 0.002
    outer = 0.006
    film = 8000.0
    ambient = 15.0
    outer_temperature = 80.0
    conduct = 78.0
    rings = []
    for radius, count in ((inner, 64), (0.0032, 48), (0.0044, 56), (outer, 72)):
        angle = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
        rings.append(np.column_stack([radius * np.cos(angle), radius * np.sin(angle)]))
    points = np.vstack(rings)
    triangulation = Delaunay(points).simplices
    centroid_radius = np.hypot(points[triangulation].mean(axis=1)[:, 0], points[triangulation].mean(axis=1)[:, 1])
    triangulation = triangulation[(centroid_radius > inner + 1e-9) & (centroid_radius < outer - 1e-9)]
    signed = _signed_area(points, triangulation)
    triangulation[signed < 0.0] = triangulation[signed < 0.0][:, ::-1]
    edge_count = {}
    for tri in triangulation:
        for i, j in ((int(tri[0]), int(tri[1])), (int(tri[1]), int(tri[2])), (int(tri[2]), int(tri[0]))):
            key = (i, j) if i < j else (j, i)
            edge_count[key] = edge_count.get(key, 0) + 1
    inner_edges = []
    outer_nodes = set()
    for i, j in (key for key, count in edge_count.items() if count == 1):
        midpoint = 0.5 * np.hypot(points[i, 0] + points[j, 0], points[i, 1] + points[j, 1])
        if midpoint < 0.5 * (inner + outer):
            inner_edges.append((i, j))
        else:
            outer_nodes.add(i)
            outer_nodes.add(j)
    mesh = {
        "points": points,
        "triangles": triangulation,
        "hole_edges": np.array(inner_edges, dtype=int),
        "flux_edges": np.zeros((0, 2), dtype=int),
        "air_edges": np.zeros((0, 2), dtype=int),
        "shank_nodes": np.array(sorted(outer_nodes), dtype=int),
    }
    spec = {
        "conductivity_ref_w_mk": conduct,
        "conductivity_softening_per_k": 0.0,
        "conductivity_ref_temperature_c": 0.0,
        "coolant_h_w_m2k": film,
        "coolant_temperature_c": ambient,
        "air_h_w_m2k": 0.0,
        "air_temperature_c": ambient,
        "shank_temperature_c": outer_temperature,
        "heat_flux_per_speed_w_m2_per_m_s": 0.0,
        "stress_free_temperature_c": ambient,
    }
    temperature = solve_temperature(mesh, spec, 1.0)
    coefficient = film * (outer_temperature - ambient) / (conduct / inner - film * np.log(inner / outer))
    offset = outer_temperature - coefficient * np.log(outer)
    radius = np.hypot(points[:, 0], points[:, 1])
    exact = offset + coefficient * np.log(np.maximum(radius, inner))
    return float(np.max(np.abs(temperature - exact)))


def _corner_index(mesh):
    points = mesh["points"]
    return int(np.argmin(points[:, 0] ** 2 + (points[:, 1] - points[:, 1].max()) ** 2))


def _sample_top(mesh, temperature, target):
    points = mesh["points"]
    for i, j in mesh["flux_edges"]:
        a = points[i]
        b = points[j]
        if abs(a[1] - target[1]) > 1e-8 or abs(b[1] - target[1]) > 1e-8:
            continue
        span = b[0] - a[0]
        if span == 0.0:
            continue
        xi = (target[0] - a[0]) / span
        if -1e-8 <= xi <= 1.0 + 1e-8:
            xi = min(1.0, max(0.0, xi))
            return float((1.0 - xi) * temperature[i] + xi * temperature[j])
    distance = np.hypot(points[:, 0] - target[0], points[:, 1] - target[1])
    return float(temperature[int(np.argmin(distance))])


def _gradients(coords):
    x, y = coords[:, 0], coords[:, 1]
    area = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))
    grads = np.array([[y[1] - y[2], y[2] - y[0], y[0] - y[1]], [x[2] - x[1], x[0] - x[2], x[1] - x[0]]]) / (2.0 * area)
    return abs(area), grads


def _conduction(coords, conduct):
    area, grads = _gradients(coords)
    return conduct * area * (grads.T @ grads), area


def _strain_matrix(grads):
    matrix = np.zeros((3, 6))
    for node in range(3):
        matrix[0, 2 * node] = grads[0, node]
        matrix[1, 2 * node + 1] = grads[1, node]
        matrix[2, 2 * node] = grads[1, node]
        matrix[2, 2 * node + 1] = grads[0, node]
    return matrix


def _plane_strain_modulus(young, poisson):
    scale = young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    return scale * np.array(
        [
            [1.0 - poisson, poisson, 0.0],
            [poisson, 1.0 - poisson, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - 2.0 * poisson)],
        ]
    )


def _scatter(rows, cols, vals, dofs, local, _stride):
    for a, dof_a in enumerate(dofs):
        for b, dof_b in enumerate(dofs):
            rows.append(int(dof_a))
            cols.append(int(dof_b))
            vals.append(float(local[a, b]))


def _edge_length(points, i, j):
    return float(np.linalg.norm(points[j] - points[i]))


def _robin(points, edges, film, ambient, rows, cols, vals, rhs):
    for i, j in edges:
        length = _edge_length(points, i, j)
        local = film * length / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
        _scatter(rows, cols, vals, (int(i), int(j)), local, 1)
        rhs[int(i)] += film * ambient * length / 2.0
        rhs[int(j)] += film * ambient * length / 2.0


def _neumann(points, edges, flux, rhs):
    for i, j in edges:
        length = _edge_length(points, i, j)
        rhs[int(i)] += flux * length / 2.0
        rhs[int(j)] += flux * length / 2.0
