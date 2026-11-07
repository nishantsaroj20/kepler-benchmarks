"""Axisymmetric shrink-fit model for a short hub on a hollow shaft.

Linear quadrilaterals, isotropic elasticity, stuck contact. The reference
solution and the truth generator each carry this file and call it on
different meshes. Tests do not import it.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def _elastic_matrix(young, poisson):
    scale = young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    normal = 1.0 - poisson
    shear = 0.5 * (1.0 - 2.0 * poisson)
    return scale * np.array(
        [
            [normal, poisson, poisson, 0.0],
            [poisson, normal, poisson, 0.0],
            [poisson, poisson, normal, 0.0],
            [0.0, 0.0, 0.0, shear],
        ]
    )


def contact_solution(problem, omega, shaft_temperature_c, sleeve_temperature_c, mesh):
    """Return nodal contact pressures and the torque capacity of one case."""
    inner = float(problem["shaft_inner_radius_m"])
    radius = float(problem["shaft_outer_radius_m"])
    outer = float(problem["sleeve_outer_radius_m"])
    shaft_length = float(problem["shaft_length_m"])
    sleeve_length = float(problem["sleeve_length_m"])
    sleeve_z0 = float(problem["sleeve_z0_m"])
    radial_interference = 0.5 * float(problem["diametral_interference_m"])
    friction = float(problem["friction_coefficient"])
    assembly_c = float(problem["assembly_temperature_c"])
    shaft = problem["shaft"]
    sleeve = problem["sleeve"]

    nr_shaft = int(mesh["nr_shaft"])
    nr_sleeve = int(mesh["nr_sleeve"])
    nz_contact = int(mesh["nz_contact"])
    nz_side = int(mesh["nz_side"])
    penalty_scale = float(mesh["penalty_scale"])

    z_below = np.linspace(0.0, sleeve_z0, nz_side + 1)
    z_mid = np.linspace(sleeve_z0, sleeve_z0 + sleeve_length, nz_contact + 1)[1:]
    z_above = np.linspace(sleeve_z0 + sleeve_length, shaft_length, nz_side + 1)[1:]
    z_nodes = np.concatenate([z_below, z_mid, z_above])
    n_axial = len(z_nodes) - 1
    contact_start = nz_side
    contact_end = nz_side + nz_contact

    r_shaft = np.linspace(inner, radius, nr_shaft + 1)
    r_sleeve = np.linspace(radius, outer, nr_sleeve + 1)
    n_shaft_r = len(r_shaft)
    n_sleeve_r = len(r_sleeve)
    n_z = len(z_nodes)
    n_shaft_nodes = n_shaft_r * n_z
    sleeve_rows = list(range(contact_start, contact_end + 1))
    sleeve_row_index = {row: index for index, row in enumerate(sleeve_rows)}
    n_nodes = n_shaft_nodes + n_sleeve_r * len(sleeve_rows)

    def shaft_node(i_r, i_z):
        return i_z * n_shaft_r + i_r

    def sleeve_node(i_r, i_z):
        return n_shaft_nodes + sleeve_row_index[i_z] * n_sleeve_r + i_r

    def coordinates(node):
        if node < n_shaft_nodes:
            i_z, i_r = divmod(node, n_shaft_r)
            return float(r_shaft[i_r]), float(z_nodes[i_z])
        local = node - n_shaft_nodes
        local_row, i_r = divmod(local, n_sleeve_r)
        return float(r_sleeve[i_r]), float(z_nodes[sleeve_rows[local_row]])

    elements = []
    for i_z in range(n_axial):
        for i_r in range(nr_shaft):
            elements.append(
                (
                    shaft_node(i_r, i_z),
                    shaft_node(i_r + 1, i_z),
                    shaft_node(i_r + 1, i_z + 1),
                    shaft_node(i_r, i_z + 1),
                    0,
                )
            )
        if contact_start <= i_z < contact_end:
            for i_r in range(nr_sleeve):
                elements.append(
                    (
                        sleeve_node(i_r, i_z),
                        sleeve_node(i_r + 1, i_z),
                        sleeve_node(i_r + 1, i_z + 1),
                        sleeve_node(i_r, i_z + 1),
                        1,
                    )
                )

    n_dof = 2 * n_nodes
    row_index = []
    col_index = []
    values = []
    load = np.zeros(n_dof)
    gauss = 1.0 / np.sqrt(3.0)
    gauss_points = (
        (-gauss, -gauss),
        (gauss, -gauss),
        (gauss, gauss),
        (-gauss, gauss),
    )
    temperature_rise = (
        shaft_temperature_c - assembly_c,
        sleeve_temperature_c - assembly_c,
    )
    materials = (shaft, sleeve)

    def add_element(dofs, stiffness, force):
        for local_a, dof_a in enumerate(dofs):
            load[dof_a] += force[local_a]
            for local_b, dof_b in enumerate(dofs):
                value = stiffness[local_a, local_b]
                if value != 0.0:
                    row_index.append(dof_a)
                    col_index.append(dof_b)
                    values.append(value)

    for n0, n1, n2, n3, material_id in elements:
        nodes = (n0, n1, n2, n3)
        xy = np.array([coordinates(node) for node in nodes])
        material = materials[material_id]
        modulus = _elastic_matrix(material["young_modulus_pa"], material["poisson_ratio"])
        thermal_strain = material["alpha_per_c"] * temperature_rise[material_id] * np.array(
            [1.0, 1.0, 1.0, 0.0]
        )
        stiffness = np.zeros((8, 8))
        force = np.zeros(8)
        for xi, eta in gauss_points:
            shape = 0.25 * np.array(
                [
                    (1.0 - xi) * (1.0 - eta),
                    (1.0 + xi) * (1.0 - eta),
                    (1.0 + xi) * (1.0 + eta),
                    (1.0 - xi) * (1.0 + eta),
                ]
            )
            shape_grad = 0.25 * np.array(
                [
                    [-(1.0 - eta), 1.0 - eta, 1.0 + eta, -(1.0 + eta)],
                    [-(1.0 - xi), -(1.0 + xi), 1.0 + xi, 1.0 - xi],
                ]
            )
            jacobian = shape_grad @ xy
            determinant = float(np.linalg.det(jacobian))
            if determinant <= 0.0:
                raise RuntimeError("element jacobian is not positive")
            physical_grad = np.linalg.solve(jacobian, shape_grad)
            radius_g = float(shape @ xy[:, 0])
            operator = np.zeros((4, 8))
            for node_i in range(4):
                operator[0, 2 * node_i] = physical_grad[0, node_i]
                operator[1, 2 * node_i + 1] = physical_grad[1, node_i]
                operator[2, 2 * node_i] = shape[node_i] / radius_g
                operator[3, 2 * node_i] = physical_grad[1, node_i]
                operator[3, 2 * node_i + 1] = physical_grad[0, node_i]
            weight = 2.0 * np.pi * radius_g * determinant
            stiffness += (operator.T @ modulus @ operator) * weight
            force += (operator.T @ modulus @ thermal_strain) * weight
            body = material["density_kg_m3"] * omega * omega * radius_g
            for node_i in range(4):
                force[2 * node_i] += shape[node_i] * body * weight
        dofs = []
        for node in nodes:
            dofs.extend((2 * node, 2 * node + 1))
        add_element(dofs, stiffness, force)

    strip = sleeve_length / nz_contact
    spring = penalty_scale * max(shaft["young_modulus_pa"], sleeve["young_modulus_pa"])
    spring *= 2.0 * np.pi * radius * strip
    for i_z in range(contact_start, contact_end + 1):
        shaft_id = shaft_node(n_shaft_r - 1, i_z)
        sleeve_id = sleeve_node(0, i_z)
        for component, jump in ((0, radial_interference), (1, 0.0)):
            shaft_dof = 2 * shaft_id + component
            sleeve_dof = 2 * sleeve_id + component
            row_index.extend((sleeve_dof, shaft_dof, sleeve_dof, shaft_dof))
            col_index.extend((sleeve_dof, shaft_dof, shaft_dof, sleeve_dof))
            values.extend((spring, spring, -spring, -spring))
            load[sleeve_dof] += spring * jump
            load[shaft_dof] -= spring * jump

    anchor = 2 * shaft_node(0, 0) + 1
    row_index.append(anchor)
    col_index.append(anchor)
    values.append(1.0e8 * spring)

    matrix = coo_matrix((values, (row_index, col_index)), shape=(n_dof, n_dof)).tocsc()
    displacement = spsolve(matrix, load)
    if not np.all(np.isfinite(displacement)):
        raise RuntimeError("displacement solution is not finite")

    pressures = []
    z_contact = []
    gaps = []
    for i_z in range(contact_start, contact_end + 1):
        shaft_id = shaft_node(n_shaft_r - 1, i_z)
        sleeve_id = sleeve_node(0, i_z)
        gap = displacement[2 * sleeve_id] - displacement[2 * shaft_id] - radial_interference
        outward_force = -spring * gap
        end_node = i_z == contact_start or i_z == contact_end
        area = 2.0 * np.pi * radius * strip * (0.5 if end_node else 1.0)
        pressures.append(outward_force / area)
        z_contact.append(float(z_nodes[i_z]))
        gaps.append(gap)

    torque = 0.0
    for i_z, pressure in zip(range(contact_start, contact_end + 1), pressures):
        end_node = i_z == contact_start or i_z == contact_end
        area = 2.0 * np.pi * radius * strip * (0.5 if end_node else 1.0)
        torque += friction * pressure * area * radius

    return {
        "pressure_pa": np.asarray(pressures, dtype=float),
        "z_m": np.asarray(z_contact, dtype=float),
        "torque_nm": float(torque),
        "gap_rms_m": float(np.sqrt(np.mean(np.square(gaps)))),
        "min_pressure_pa": float(np.min(pressures)),
    }


def station_pressures(problem, solved):
    z0 = float(problem["sleeve_z0_m"])
    length = float(problem["sleeve_length_m"])
    fractions = (0.1, 0.3, 0.5, 0.7, 0.9)
    return [
        float(np.interp(z0 + fraction * length, solved["z_m"], solved["pressure_pa"]))
        for fraction in fractions
    ]


def limit_speed(problem, mesh):
    """Bisect the speed whose hot torque capacity equals the required torque."""
    required = float(problem["required_torque_nm"])
    low = 1000.0
    high = 3000.0
    shaft_c = float(problem["operating_shaft_temperature_c"])
    sleeve_c = float(problem["operating_sleeve_temperature_c"])
    low_case = contact_solution(problem, low, shaft_c, sleeve_c, mesh)
    high_case = contact_solution(problem, high, shaft_c, sleeve_c, mesh)
    if low_case["min_pressure_pa"] <= 0.0 or high_case["min_pressure_pa"] <= 0.0:
        raise RuntimeError("contact unloads inside the speed bracket")
    if not (low_case["torque_nm"] > required > high_case["torque_nm"]):
        raise RuntimeError("required torque is not bracketed by 1000 and 3000 rad/s")
    for _ in range(30):
        mid = 0.5 * (low + high)
        torque = contact_solution(problem, mid, shaft_c, sleeve_c, mesh)["torque_nm"]
        if torque > required:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def predict(problem, mesh):
    assembly_c = float(problem["assembly_temperature_c"])
    assembly = contact_solution(problem, 0.0, assembly_c, assembly_c, mesh)
    if assembly["min_pressure_pa"] <= 0.0:
        raise RuntimeError("assembly contact is not compressive")
    return {
        "assembly_pressure_pa": station_pressures(problem, assembly),
        "assembly_torque_capacity_nm": assembly["torque_nm"],
        "limit_speed_rad_s": limit_speed(problem, mesh),
    }
