"""
Coordinate convention:
    X = forward
    Y = left / outboard
    Z = upward

All coordinates and lengths are in millimetres.
"""

import numpy as np
import matplotlib.pyplot as plt
import os


# ============================================================
# 1. USER INPUTS
# ============================================================

# Coordinate helper:
# Each point is entered as [X, Y, Z] in millimetres.
# Coordinates are for the left side of the vehicle. The right side is mirrored automatically.

hardpoints = {
    # Lower control arm chassis pickups
    "LCA_front_inner": np.array([200.0, 18.0, 250.0]),
    "LCA_rear_inner": np.array([-400.0, 15.0, 259.0]),

    # Upper control arm chassis pickups
    "UCA_front_inner": np.array([28.2, 162.0, 442.6]),
    "UCA_rear_inner": np.array([-340.1, 196.2, 437.5]),

    # Ball-joint positions at the design position
    "lower_ball_joint": np.array([-2.5, 700.5, 267.3]),
    "upper_ball_joint": np.array([-24, 654.8, 474.7]),

    # Wheel Center in upright
    "wheel_center": np.array([-20.0, 788.2, 350.0,]),

    # Wheel Hub axis point
    "hub_axis_point": np.array([-20.0, 800.0, 350.0]),

    # Tie rods
    "tie_rod_inner": np.array([149.4, 170.0, 446.7]),
    "tie_rod_outer": np.array([73.8, 683.8, 458.3]),

    # Wheel radius
    "wheel_radius_mm": float(350.0),

}

# LCA angle sweep settings, in degrees
ANGLE_START_DEG = -2.5
ANGLE_END_DEG = 2.5
NUMBER_OF_POSITIONS = 41

# Angle used for the 3D suspension picture
DISPLAY_ANGLE_DEG = 0.0

# Figure export settings
FIGURE_DPI = 300
OUTPUT_FOLDER = "portfolio_figures"

# ============================================================
# 2. BASIC VECTOR FUNCTIONS
# ============================================================

def vector_length(vector):

    return np.linalg.norm(vector)


def unit_vector(vector):

    length = vector_length(vector)

    if length < 1e-12:
        raise ValueError("Cannot normalize a vector with zero length.")

    return vector / length


def distance(point_a, point_b):

    return vector_length(point_b - point_a)


# ============================================================
# 3. MIRROR FUNCTION
# ============================================================

def mirror_point(point):

    mirrored = point.copy()
    mirrored[1] = -1.0 * mirrored[1]

    return mirrored

# ============================================================
# 4. MIRROR ALL HARDPOINTS
# ============================================================

def mirror_hardpoints(left_hardpoints):
    
    mirrored = {}

    for key, value in left_hardpoints.items():

        if isinstance(value, np.ndarray):
            mirrored[key] = mirror_point(value)

        else:
            mirrored[key] = value

    return mirrored

# ============================================================
# 5. ROTATION ABOUT AN ARBITRARY AXIS
# ============================================================

def rotate_point_about_axis(point, axis_point_1, axis_point_2, angle_rad):
    """
    Parameters
    ----------
    point:
        Point to rotate.

    axis_point_1:
        First point defining the rotation axis.

    axis_point_2:
        Second point defining the rotation axis.

    angle_rad:
        Rotation angle in radians.

    Returns
    -------
    rotated_point:
        New point position after rotation.

    Notes
    -----
    This function uses Rodrigues' rotation formula.
    """

    # Unit vector along the rotation axis
    axis_direction = unit_vector(axis_point_2 - axis_point_1) 

    # Move the coordinate system so that axis_point_1 is the origin
    relative_point = point - axis_point_1

    # Rodrigues' rotation formula
    rotated_relative_point = (
        relative_point * np.cos(angle_rad)
        + np.cross(axis_direction, relative_point) * np.sin(angle_rad)
        + axis_direction
        * np.dot(axis_direction, relative_point)
        * (1.0 - np.cos(angle_rad))
    )

    # Move back to the original global coordinate system
    rotated_point = axis_point_1 + rotated_relative_point

    return rotated_point

# ============================================================
# 6. INTERSECTION OF THREE SPHERES
# ============================================================

def three_sphere_intersection(
    center_1,
    radius_1,
    center_2,
    radius_2,
    center_3,
    radius_3,
):
    """
    Find the intersection points of three spheres.

    Each sphere is defined by:

        ||P - center|| = radius

    In this suspension model, the unknown point P is the upper
    ball joint.

    Sphere 1:
        Center = UCA front inner pickup
        Radius = front leg length of upper wishbone

    Sphere 2:
        Center = UCA rear inner pickup
        Radius = rear leg length of upper wishbone

    Sphere 3:
        Center = current lower ball joint
        Radius = upright length

    Returns
    -------
    solution_1, solution_2:
        The two possible intersection points.

    Raises
    ------
    ValueError:
        If the three spheres have no valid intersection.
    """

    # Unit direction from center 1 to center 2
    ex_vector = center_2 - center_1
    d = vector_length(ex_vector)

    if d < 1e-12:
        raise ValueError(
            "Sphere centers 1 and 2 are at the same location."
        )

    ex = ex_vector / d

    # Vector from center 1 to center 3
    center_1_to_3 = center_3 - center_1

    # Projection of center 1 -> center 3 onto ex
    i = np.dot(ex, center_1_to_3)

    # Construct a second perpendicular basis direction
    temporary_vector = center_1_to_3 - i * ex
    temporary_length = vector_length(temporary_vector)

    if temporary_length < 1e-12:
        raise ValueError(
            "The three sphere centers are collinear or nearly collinear."
        )

    ey = temporary_vector / temporary_length

    # Third mutually perpendicular direction
    ez = np.cross(ex, ey)

    # Coordinates of center 3 in this local coordinate system
    j = np.dot(ey, center_1_to_3)

    # Trilateration equations
    x = (
        radius_1**2
        - radius_2**2
        + d**2
    ) / (2.0 * d)

    y = (
        radius_1**2
        - radius_3**2
        + i**2
        + j**2
        - 2.0 * i * x
    ) / (2.0 * j)

    z_squared = radius_1**2 - x**2 - y**2

    # Small negative values can occur because of floating-point rounding
    if z_squared < -1e-8:
        raise ValueError(
            "The three spheres do not intersect. "
            f"Calculated z^2 = {z_squared:.6f}"
        )

    z = np.sqrt(max(z_squared, 0.0))

    # Base point in the plane of the three sphere centers
    base_point = center_1 + x * ex + y * ey

    # Two possible solutions, one on either side of the plane
    solution_1 = base_point + z * ez
    solution_2 = base_point - z * ez

    return solution_1, solution_2

# ============================================================
# 7. INTERSECT LINE WITH HORIZONTAL PLANE
# ============================================================

def intersect_line_with_horizontal_plane(line_point, line_direction, plane_z):

    dz = line_direction[2]
    if abs(dz) < 1e-12:
        raise ValueError(
            "Cannot intersect line with ground plane because "
            "the line is parallel or nearly parallel to the plane."
        )

    parameter = (plane_z - line_point[2]) / dz
    intersection = line_point + parameter*line_direction

    return intersection

# ============================================================
# 8. INTERSECT TWO LINES
# ============================================================
def intersect_lines_in_2d(line_point_1, line_direction_1, line_point_2, line_direction_2):

    A = np.array([
        [line_direction_1[0], -line_direction_2[0]], 
        [line_direction_1[1], -line_direction_2[1]]
        ])

    b = np.array([
        [line_point_2[0]-line_point_1[0]], 
        [line_point_2[1]-line_point_1[1]]
        ])

    determinant = np.linalg.det(A)
    if abs(determinant) < 1e-12:
        raise ValueError("The two lines are parallel or nearly parallel.")

    # Find parameter for one of the line
    x = np.linalg.solve(A, b)
    t = x[0]

    intersection = line_point_1 + t*line_direction_1

    return intersection

# ============================================================
# 9. CONVERT 3D POINT TO Y-Z AXES
# ============================================================

def point_to_yz(point_3d):

    return np.array([
        point_3d[1],
        point_3d[2],
    ])

# ============================================================
# 10. FIND PERPENDICULAR DIRECTION OF 2D VECTOR
# ============================================================

def perpendicular_direction_2d(vector_2d):

    if vector_length(vector_2d) < 1e-12:
        raise ValueError(
            "Cannot calculate a perpendicular direction from "
            "a zero-length 2D vector."
        )

    return np.array([
        -vector_2d[1],
        vector_2d[0],
    ])

# ============================================================
# 11. CALCULATE LOWEST WHEEL POINT
# ============================================================

def calculate_lowest_wheel_point(wheel_center, hub_axis_point, wheel_radius):

    wheel_axis = unit_vector(hub_axis_point - wheel_center)
    global_down = np.array([0.0,0.0,-1.0])
    wheel_plane_vector = unit_vector(global_down - np.dot(global_down, wheel_axis)*wheel_axis)

    lowest_wheel_point = wheel_center + wheel_radius*wheel_plane_vector

    return lowest_wheel_point

# ============================================================
# 12. CREATE THE SUSPENSION MODEL
# ============================================================

def create_suspension_model(input_hardpoints):

    # Copy the arrays so that the user's original dictionary is not modified accidentally.
    model = {
        name: (
            value.astype(float).copy()
            if isinstance(value, np.ndarray)
            else float(value)
        )
        for name, value in input_hardpoints.items()
    }

    # Determine suspension side
    model["side_sign"] = float(np.sign(model["wheel_center"][1]))
    
    if model["side_sign"] == 0.0:
        raise ValueError(
            "Cannot determine suspension side because "
            "wheel-center Y coordinate is zero."
        )                    

    # Calculate fixed link lengths from the design position
    model["LCA_front_length"] = distance(
        model["LCA_front_inner"],
        model["lower_ball_joint"],
    )

    model["LCA_rear_length"] = distance(
        model["LCA_rear_inner"],
        model["lower_ball_joint"],
    )

    model["UCA_front_length"] = distance(
        model["UCA_front_inner"],
        model["upper_ball_joint"],
    )

    model["UCA_rear_length"] = distance(
        model["UCA_rear_inner"],
        model["upper_ball_joint"],
    )

    model["upright_length"] = distance(
        model["lower_ball_joint"],
        model["upper_ball_joint"],
    )

    model["tie_rod_length"] = distance(
        model["tie_rod_inner"],
        model["tie_rod_outer"],
    )

    model["lbj_to_tro_length"] = distance(
        model["lower_ball_joint"],
        model["tie_rod_outer"],
    )

    model["ubj_to_tro_length"] = distance(
        model["upper_ball_joint"],
        model["tie_rod_outer"],
    )

    model["wc_to_lbj_length"] = distance(
        model["wheel_center"],
        model["lower_ball_joint"],
    )

    model["wc_to_hub_length"] = distance(
        model["wheel_center"],
        model["hub_axis_point"],
    )

    # Build local upright axes relative to LBJ
    lbj = model["lower_ball_joint"]
    ubj = model["upper_ball_joint"]
    wheel_center = model["wheel_center"]
    tie_rod_outer = model["tie_rod_outer"]
    hub_axis_point = model["hub_axis_point"]

    upright_z = unit_vector(ubj - lbj) #local upright z coordinate

    tie_vector = tie_rod_outer-lbj
    tie_vector_projection = tie_vector-(np.dot(tie_vector, upright_z))*upright_z

    upright_y = unit_vector(tie_vector_projection) #local upright y coordinate

    upright_x = unit_vector(np.cross(upright_y, upright_z)) #local upright x coordinate

    # Local coordinates of Wheel Center
    wc_relative = wheel_center-lbj

    local_x_wc = np.dot(wc_relative, upright_x)
    local_y_wc = np.dot(wc_relative, upright_y)
    local_z_wc = np.dot(wc_relative, upright_z)

    model["wheel_center_local"] = np.array([local_x_wc, local_y_wc, local_z_wc])

    # Local coordinates of hub axis point
    hap_relative = hub_axis_point - lbj

    local_x_hap = np.dot(hap_relative, upright_x)
    local_y_hap = np.dot(hap_relative, upright_y)
    local_z_hap = np.dot(hap_relative, upright_z)

    model["hub_axis_point_local"] = np.array([local_x_hap, local_y_hap, local_z_hap])        

    return model

# ============================================================
# 13. SOLVE ONE SUSPENSION POSITION
# ============================================================

def solve_suspension_position(model, lca_angle_deg, ubj_reference, tie_rod_reference):
    """
    Solve the suspension position for one LCA rotation angle.

    Procedure
    ---------
    1. Rotate the lower ball joint around the LCA inboard axis.
    2. Use the UCA lengths and upright length to locate the upper ball joint.
    3. Select the sphere-intersection solution closest to the design upper-ball-joint position.
    4. Solve the position of the outer tie rod using the same method.
    4. Use the new ball-joint and outer tie rod positions to construct a new upright frame to find the coordinate of the new wheel center.
    """

    # ------LBJ & UBJ------
    angle_rad = np.deg2rad(lca_angle_deg)

    # Rotate the lower ball joint around the lower-arm pivot axis
    current_lower_ball_joint = rotate_point_about_axis(
        point=model["lower_ball_joint"],
        axis_point_1=model["LCA_front_inner"],
        axis_point_2=model["LCA_rear_inner"],
        angle_rad=angle_rad,
    )

    # Calculate the two mathematically possible locations of the upper ball joint
    upper_solution_1, upper_solution_2 = three_sphere_intersection(
        center_1=model["UCA_front_inner"],
        radius_1=model["UCA_front_length"],

        center_2=model["UCA_rear_inner"],
        radius_2=model["UCA_rear_length"],

        center_3=current_lower_ball_joint,
        radius_3=model["upright_length"],
    )

    # Select the solution closest to the reference position.
    # This prevents selecting the mirrored assembly configuration.
    distance_to_solution_1 = distance(
        upper_solution_1,
        ubj_reference,
    )

    distance_to_solution_2 = distance(
        upper_solution_2,
        ubj_reference,
    )

    if distance_to_solution_1 <= distance_to_solution_2:
        current_upper_ball_joint = upper_solution_1
    else:
        current_upper_ball_joint = upper_solution_2

    # ------TIE ROD------
    tie_rod_1, tie_rod_2 = three_sphere_intersection(
            center_1=model["tie_rod_inner"],
            radius_1=model["tie_rod_length"],
    
            center_2=current_lower_ball_joint,
            radius_2=model["lbj_to_tro_length"],
    
            center_3=current_upper_ball_joint,
            radius_3=model["ubj_to_tro_length"],
    )

    # Select the solution closest to the reference position.
    distance_to_solution_1 = distance(
        tie_rod_1,
        tie_rod_reference,
    )
    
    distance_to_solution_2 = distance(
        tie_rod_2,
        tie_rod_reference,
    )

    if distance_to_solution_1 <= distance_to_solution_2:
            current_tie_rod = tie_rod_1
    else:
        current_tie_rod = tie_rod_2

    # ------WHEEL CENTER------
    lbj = current_lower_ball_joint
    ubj = current_upper_ball_joint
    tie_rod_outer = current_tie_rod
    wc_local = model["wheel_center_local"]
    hub_local = model["hub_axis_point_local"]

    # Construct new upright coordinates
    new_upright_z = unit_vector(ubj - lbj) #local upright z coordinate
    
    tie_vector = tie_rod_outer-lbj
    tie_vector_projection = tie_vector-(np.dot(tie_vector, new_upright_z))*new_upright_z
    
    new_upright_y = unit_vector(tie_vector_projection) #local upright y coordinate
    
    new_upright_x = unit_vector(np.cross(new_upright_y, new_upright_z)) #local upright x coordinate

    # Find current wheel center
    current_wheel_center = lbj + wc_local[0]*new_upright_x + wc_local[1]*new_upright_y + wc_local[2]*new_upright_z

    # ------ HUB AXIS POINT ------
    current_hub_axis_point = lbj + hub_local[0]*new_upright_x + hub_local[1]*new_upright_y + hub_local[2]*new_upright_z

    # Construct results dictionary
    results = {
        "lca_angle_deg": lca_angle_deg,
        "lower_ball_joint": current_lower_ball_joint,
        "upper_ball_joint": current_upper_ball_joint,
        "tie_rod_outer": current_tie_rod,
        "wheel_center": current_wheel_center,
        "hub_axis_point": current_hub_axis_point,
    }

    return results

# ============================================================
# 14. CASTER AND KPI
# ============================================================

def calculate_alignment_angles(model, results):

    lbj = results["lower_ball_joint"]
    ubj = results["upper_ball_joint"]

    side_sign = model["side_sign"]

    steering_axis = ubj - lbj

    dx = steering_axis[0]
    dy = steering_axis[1]
    dz = steering_axis[2]

    # X and Z do not change under left-right reflection.
    caster_deg = np.degrees(
        np.arctan2(-dx, dz)
    )

    # Positive KPI when top of the steering axis points inboard
    kpi_deg = np.degrees(
        np.arctan2(-side_sign * dy, dz)
    )

    return {
        "caster_deg": caster_deg,
        "kpi_deg": kpi_deg,
    }

# ============================================================
# 15. WHEEL TRAVEL
# ============================================================

def calculate_wheel_center_travel(model, results):

    wheel_center_travel = results["wheel_center"] - model["wheel_center"]

    wc_travel_x = wheel_center_travel[0]
    wc_travel_y = wheel_center_travel[1]
    wc_travel_z = wheel_center_travel[2]

    wheel_travel = {
        "wc_travel_x": wc_travel_x,
        "wc_travel_y": wc_travel_y,
        "wc_travel_z": wc_travel_z,
    }

    return wheel_travel

# ============================================================
# 16. CAMBER & TOE
# ============================================================

def calculate_camber_toe(model, results):

    wheel_center = results["wheel_center"]
    hub_axis_point = results["hub_axis_point"]

    side_sign = model["side_sign"]

    wheel_axis = hub_axis_point - wheel_center

    hx = wheel_axis[0]
    hy = wheel_axis[1]
    hz = wheel_axis[2]

    # Front-view projection onto the Y-Z plane
    # Positive camber is when top of the tire points outboard
    camber_deg = np.degrees(np.arctan2(-hz, side_sign * hy))

    # Top-view projection onto the X-Y plane
    # Positive toe means TOE IN
    toe_deg = np.degrees(np.arctan2(hx, side_sign * hy))

    return {
        "camber_deg": camber_deg,
        "toe_deg": toe_deg,
    }   

# ============================================================
# 17. CALCULATE SCRUB RADIUS AND MECHANICAL TRAIL
# ============================================================

def calculate_sr_mt(model, results):

    ubj = results["upper_ball_joint"]   
    lbj = results["lower_ball_joint"]
    wheel_center = results["wheel_center"]
    hub_axis_point = results["hub_axis_point"]
    wheel_radius = model["wheel_radius_mm"]

    side_sign = model["side_sign"]

    lowest_wheel_point = calculate_lowest_wheel_point(wheel_center, hub_axis_point, wheel_radius)
    ground_z = lowest_wheel_point[2]
    steering_axis = unit_vector(lbj-ubj)

    # Set coordinates of steering axis projection and contact patch center
    steering_axis_ground = intersect_line_with_horizontal_plane(ubj, steering_axis, ground_z)
    contact_patch_center = lowest_wheel_point

    # Calculate SCRUB RADIUS
    # SIGN CONVENTION: 
    # POSITIVE scrub radius when contact patch center is OUTBOARD of projected steering axis
    scrub_radius_mm = side_sign * (contact_patch_center[1] - steering_axis_ground[1])

    # Calculate MECHANICAL TRAIL   
    # SIGN CONVENTION:
    # POSITIVE mechanical trail when projected steering axis is AHEAD of contact patch center
    mechanical_trail_mm = steering_axis_ground[0] - contact_patch_center[0]

    return {
        "ground_z": ground_z,
        "steering_axis_ground": steering_axis_ground,
        "contact_patch_center": contact_patch_center,
        "scrub_radius_mm": scrub_radius_mm,
        "mechanical_trail_mm": mechanical_trail_mm,
    }

# ============================================================
# 18. FIND FRONT VIEW INSTANT CENTER   
# ============================================================
def calculate_front_view_instant_center(model, current_results, angle_step_deg=0.001):

    if angle_step_deg <= 0.0:
        raise ValueError("angle_step_deg must be greater than zero.")

    current_angle_deg = current_results["lca_angle_deg"]
    current_ubj = current_results["upper_ball_joint"]
    current_tro = current_results["tie_rod_outer"]

    # Solve positions slightly above and below the current LCA angle.
    results_plus = solve_suspension_position(
        model=model,
        lca_angle_deg=current_angle_deg + angle_step_deg,
        ubj_reference=current_ubj,
        tie_rod_reference=current_tro,
    )

    results_minus = solve_suspension_position(
        model=model,
        lca_angle_deg=current_angle_deg - angle_step_deg,
        ubj_reference=current_ubj,
        tie_rod_reference=current_tro,
    )

    # Current front-view joint coordinates
    current_lbj_yz = point_to_yz(current_results["lower_ball_joint"])

    current_ubj_yz = point_to_yz(current_results["upper_ball_joint"])

    # Approximate front-view joint velocity directions using central finite differences.
    lbj_displacement_yz = (
        point_to_yz(results_plus["lower_ball_joint"])
        - point_to_yz(results_minus["lower_ball_joint"])
    )

    ubj_displacement_yz = (
        point_to_yz(results_plus["upper_ball_joint"])
        - point_to_yz(results_minus["upper_ball_joint"])
    )

    if vector_length(lbj_displacement_yz) < 1e-12:
        raise ValueError("The projected LBJ displacement is too small to calculate an instant center.")

    if vector_length(ubj_displacement_yz) < 1e-12:
        raise ValueError("The projected UBJ displacement is too small to calculate an instant center.")

    # The line from each joint toward the instant center is perpendicular to that joint's instantaneous velocity.
    lbj_normal_direction = perpendicular_direction_2d(lbj_displacement_yz)

    ubj_normal_direction = perpendicular_direction_2d(ubj_displacement_yz)

    # Intersect the two instantaneous-normal lines.
    
    instant_center_yz = intersect_lines_in_2d(
        line_point_1=current_lbj_yz,
        line_direction_1=lbj_normal_direction,
        line_point_2=current_ubj_yz,
        line_direction_2=ubj_normal_direction,
    )

    construction = {
        "instant_center_yz": instant_center_yz,
        "lbj_yz": current_lbj_yz,
        "ubj_yz": current_ubj_yz,
        "lbj_displacement_yz": lbj_displacement_yz,
        "ubj_displacement_yz": ubj_displacement_yz,
        "lbj_normal_direction": lbj_normal_direction,
        "ubj_normal_direction": ubj_normal_direction,
        "angle_step_deg": angle_step_deg,
    }

    return instant_center_yz, construction

# ============================================================
# 19. FIND FRONT AXLE ROLL CENTER
# ============================================================

def calculate_roll_center(left_model, left_results, right_model, right_results):

    # Find front-view instant centers of both sides of the suspension
    left_ic_yz, left_ic_construction = calculate_front_view_instant_center(left_model, left_results)
    right_ic_yz, right_ic_construction = calculate_front_view_instant_center(right_model, right_results)

    # Find front_view contact patch centers of both wheels
    left_contact_patch_yz = point_to_yz(left_results["contact_patch_center"])
    right_contact_patch_yz = point_to_yz(right_results["contact_patch_center"])

    # Find vector from ic to cpc
    left_force_line_direction = left_contact_patch_yz - left_ic_yz
    right_force_line_direction = right_contact_patch_yz - right_ic_yz

    roll_center_yz = intersect_lines_in_2d(left_ic_yz, left_force_line_direction, right_ic_yz, right_force_line_direction)    

    construction = {
        "roll_center_yz": roll_center_yz,
        "left_instant_center_yz": left_ic_yz,
        "right_instant_center_yz": right_ic_yz,
        "left_contact_patch_yz": left_contact_patch_yz,
        "right_contact_patch_yz": right_contact_patch_yz,
        "left_force_line_direction": left_force_line_direction,
        "right_force_line_direction": right_force_line_direction,
        "left_ic_construction": left_ic_construction,
        "right_ic_construction": right_ic_construction,
        }

    return roll_center_yz, construction

# ============================================================
# 20. CALCULATE ROLL CENTER SWEEP
# ============================================================
def calculate_roll_center_sweep(left_sweep_results, right_sweep_results, left_model, right_model):

    if len(left_sweep_results) != len(right_sweep_results):
        raise ValueError(
            "Left and right sweep-result lists have "
            "different lengths."
        )

    roll_center_sweep = []

    for left_result, right_result in zip(left_sweep_results, right_sweep_results):

        # Calculate the front roll center.
        roll_center_yz, construction = calculate_roll_center(left_model, left_result, right_model, right_result)

        # Average vertical wheel-center travel.
        average_wheel_travel_z = 0.5 * (left_result["wheel_travel_z_mm"] + right_result["wheel_travel_z_mm"])

        # Use one common virtual ground plane.
        left_ground_z = left_result["ground_z"]
        right_ground_z = right_result["ground_z"]

        common_ground_z = 0.5 * (left_ground_z+ right_ground_z)

        # The roll-center Z coordinate is global.
        roll_center_global_z = roll_center_yz[1]

        # This is the more useful ground-relative value.
        roll_center_height = (roll_center_global_z - common_ground_z)

        # Find roll center y
        if roll_center_yz[0] > 1e-8:
            roll_center_y = roll_center_yz[0]
        else:
            roll_center_y = np.array([0])

        roll_center_sweep.append({
            "average_wheel_travel_z_mm":
                average_wheel_travel_z,

            "left_wheel_travel_z_mm":
                left_result["wheel_travel_z_mm"],

            "right_wheel_travel_z_mm":
                right_result["wheel_travel_z_mm"],

            "roll_center_y_mm":
                roll_center_y,

            "roll_center_global_z_mm":
                roll_center_global_z,

            "common_ground_z_mm":
                common_ground_z,

            "roll_center_height_mm":
                roll_center_height,

            "left_lca_angle_deg":
                left_result["lca_angle_deg"],

            "right_lca_angle_deg":
                right_result["lca_angle_deg"],

            "construction":
                construction,
        })

    roll_center_sweep.sort(
        key=lambda result:
            result["average_wheel_travel_z_mm"]
    )

    return roll_center_sweep

# ============================================================
# 21. ADD DERIVED GEOMETRY
# ============================================================

def add_derived_geometry(model, position):

    alignment = calculate_alignment_angles(model, position)
    wheel_travel = calculate_wheel_center_travel(model,position)
    camber_toe = calculate_camber_toe(model, position)
    sr_mt = calculate_sr_mt(model,position)

    # Add alignment values
    position["caster_deg"] = alignment["caster_deg"]
    position["kpi_deg"] = alignment["kpi_deg"]

    # Add wheel-centre travel
    position["wheel_travel_x_mm"] = (wheel_travel["wc_travel_x"])
    position["wheel_travel_y_mm"] = (wheel_travel["wc_travel_y"])
    position["wheel_travel_z_mm"] = (wheel_travel["wc_travel_z"])

    # Add wheel angles
    position["camber_deg"] = camber_toe["camber_deg"]
    position["toe_deg"] = camber_toe["toe_deg"]

    # Add steering-ground geometry
    position["scrub_radius_mm"] = (sr_mt["scrub_radius_mm"])
    position["mechanical_trail_mm"] = (sr_mt["mechanical_trail_mm"])

    # Add useful construction coordinates
    position["contact_patch_center"] = (sr_mt["contact_patch_center"])
    position["steering_axis_ground"] = (sr_mt["steering_axis_ground"])
    position["ground_z"] = (sr_mt["ground_z"])

    return position

# ============================================================
# 22. CHECK THE RIGID-LINK CONSTRAINTS
# ============================================================

def calculate_constraint_errors(model, results):
    """
    Calculate the change in each rigid-link length.

    Ideally, every error should be approximately zero.
    """

    lower_ball_joint = results["lower_ball_joint"]
    upper_ball_joint = results["upper_ball_joint"]
    tie_rod_outer = results["tie_rod_outer"]
    wheel_center = results["wheel_center"]
    hub_axis_point = results["hub_axis_point"]

    errors = {
        "LCA_front_error": (
            distance(model["LCA_front_inner"], lower_ball_joint)
            - model["LCA_front_length"]
        ),

        "LCA_rear_error": (
            distance(model["LCA_rear_inner"], lower_ball_joint)
            - model["LCA_rear_length"]
        ),

        "UCA_front_error": (
            distance(model["UCA_front_inner"], upper_ball_joint)
            - model["UCA_front_length"]
        ),

        "UCA_rear_error": (
            distance(model["UCA_rear_inner"], upper_ball_joint)
            - model["UCA_rear_length"]
        ),

        "tie_rod_error": (
            distance(model["tie_rod_inner"], tie_rod_outer)
            - model["tie_rod_length"]
        ),

        "upright_error": (
            distance(lower_ball_joint, upper_ball_joint)
            - model["upright_length"]
        ),

        "LBJ to WC error": (
            distance(lower_ball_joint, wheel_center)
            - model["wc_to_lbj_length"]
        ),

        "WC to hub error": (
            distance(wheel_center, hub_axis_point)
            - model["wc_to_hub_length"]
        ),
    }

    return errors

# ============================================================
# 23. VERIFY SYMMETRY
# ============================================================

def check_mirrored_positions(left_results, right_results):
    """
    Check whether two solved positions are reflections across Y=0.
    """

    point_names = [
        "lower_ball_joint",
        "upper_ball_joint",
        "tie_rod_outer",
        "wheel_center",
        "hub_axis_point",
    ]

    maximum_error = 0.0

    print("\nMirrored coordinate check")
    print("-------------------------------------------")

    for point_name in point_names:
        left_point = left_results[point_name]
        right_point = right_results[point_name]

        expected_right = mirror_point(left_point)
        error_vector = right_point - expected_right

        point_error = np.linalg.norm(error_vector)

        maximum_error = max(
            maximum_error,
            point_error,
        )

        print(
            f"{point_name:22s}: "
            f"{point_error:.9e} mm"
        )

    print(
        f"Maximum mirror error: "
        f"{maximum_error:.9e} mm"
    )

    return maximum_error

# ============================================================
# 24. RUN A COMPLETE ANGLE SWEEP
# ============================================================

def run_angle_sweep(model, angles_deg):
    """
    Solve the suspension at every requested LCA angle.

    The sweep starts from the reference position and proceeds separately toward positive and negative LCA rotations.
    This helps preserve the correct sphere-intersection branch.
    """

    positive_angles = np.sort(
        angles_deg[angles_deg >= 0.0]
    )

    negative_angles = np.sort(
        angles_deg[angles_deg < 0.0]
    )[::-1]

    solved_positions = []

    # Sweep from zero toward positive rotation
    previous_upper_ball_joint = (model["upper_ball_joint"])
    previous_tro = (model["tie_rod_outer"])

    for angle_deg in positive_angles:
        try:
            position = solve_suspension_position(
                model,
                angle_deg,
                previous_upper_ball_joint,
                previous_tro
            )

            # Solve for derived geometries
            position = add_derived_geometry(model,position)

            solved_positions.append(position)

            previous_upper_ball_joint = (position["upper_ball_joint"])
            previous_tro = (position["tie_rod_outer"])

        except ValueError as error:
            print(
                f"Could not solve LCA angle "
                f"{angle_deg:.2f} deg: {error}"
            )
            break

    # Restart from the reference geometry and sweep negative
    previous_upper_ball_joint = (model["upper_ball_joint"])
    previous_tro = (model["tie_rod_outer"])

    for angle_deg in negative_angles:
        try:
            position = solve_suspension_position(
                model,
                angle_deg,
                previous_upper_ball_joint,
                previous_tro
            )

            # Solve for derived geometries
            position = add_derived_geometry(model,position)

            solved_positions.append(position)

            previous_upper_ball_joint = (position["upper_ball_joint"])
            previous_tro = (position["tie_rod_outer"])

        except ValueError as error:
            print(
                f"Could not solve LCA angle "
                f"{angle_deg:.2f} deg: {error}"
            )
            break

    solved_positions.sort(
        key= lambda position
        :position["lca_angle_deg"]
    )

    return solved_positions

# ============================================================
# 25. PLOT SUSPENSION KINEMATIC CURVES
# ============================================================

def plot_kinematic_curves(model, sweep_results):

    if len(sweep_results) == 0:
        raise ValueError(
            "Cannot plot kinematic curves because sweep_results "
            "is empty."
        )

    # --------------------------------------------------------
    # Extract sweep values
    # --------------------------------------------------------

    wheel_travel_z = np.array([
        result["wheel_travel_z_mm"]
        for result in sweep_results
    ])

    camber_deg = np.array([
        result["camber_deg"]
        for result in sweep_results
    ])

    toe_deg = np.array([
        result["toe_deg"]
        for result in sweep_results
    ])

    caster_deg = np.array([
        result["caster_deg"]
        for result in sweep_results
    ])

    kpi_deg = np.array([
        result["kpi_deg"]
        for result in sweep_results
    ])

    scrub_radius_mm = np.array([
        result["scrub_radius_mm"]
        for result in sweep_results
    ])

    mechanical_trail_mm = np.array([
        result["mechanical_trail_mm"]
        for result in sweep_results
    ])

    # --------------------------------------------------------
    # Sort all data by wheel travel
    # --------------------------------------------------------

    sort_order = np.argsort(wheel_travel_z)

    wheel_travel_z = wheel_travel_z[sort_order]
    camber_deg = camber_deg[sort_order]
    toe_deg = toe_deg[sort_order]
    caster_deg = caster_deg[sort_order]
    kpi_deg = kpi_deg[sort_order]
    scrub_radius_mm = scrub_radius_mm[sort_order]
    mechanical_trail_mm = mechanical_trail_mm[sort_order]

    # --------------------------------------------------------
    # Identify the design position at zero wheel travel and
    # highlight its kinematic values on each subplot
    # --------------------------------------------------------

    design_index = int(np.argmin(np.abs(wheel_travel_z)))
    design_wheel_travel_z = wheel_travel_z[design_index]

    design_values = {
        "camber_deg": camber_deg[design_index],
        "toe_deg": toe_deg[design_index],
        "caster_deg": caster_deg[design_index],
        "kpi_deg": kpi_deg[design_index],
        "scrub_radius_mm": scrub_radius_mm[design_index],
        "mechanical_trail_mm": mechanical_trail_mm[design_index],
    }

    def annotate_design_value(ax, x_value, y_value, text_value):

        x_offset = 12 if x_value >= 0 else -12
        y_offset = 10 if y_value >= ax.get_ylim()[0] + 0.5 * (ax.get_ylim()[1] - ax.get_ylim()[0]) else -10

        ax.annotate(
            text_value,
            (x_value, y_value),
            xytext=(x_offset, y_offset),
            textcoords="offset points",
            color="red",
            fontsize=9,
            weight="bold",
            ha="left" if x_value >= 0 else "right",
            va="bottom" if y_offset > 0 else "top",
            bbox=dict(
                boxstyle="round,pad=0.25",
                facecolor="white",
                edgecolor="red",
                alpha=0.85,
            ),
            arrowprops=dict(
                arrowstyle="-",
                color="red",
                lw=1,
            ),
            zorder=6,
        )

    # --------------------------------------------------------
    # Create the 3-row by 2-column figure
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        nrows=3,
        ncols=2,
        figsize=(13, 12),
        sharex=True,
    )

    # --------------------------------------------------------
    # Row 1: Camber and toe
    # --------------------------------------------------------

    axes[0, 0].plot(
        wheel_travel_z,
        camber_deg,
        marker="o",
        markersize=4,
    )

    axes[0, 0].scatter(
        [design_wheel_travel_z],
        [design_values["camber_deg"]],
        color="red",
        s=70,
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
    )
    annotate_design_value(
        axes[0, 0],
        design_wheel_travel_z,
        design_values["camber_deg"],
        f"{design_values['camber_deg']:.2f}°",
    )

    axes[0, 0].set_title("Camber vs Wheel Travel")

    axes[0, 0].set_ylabel("Camber (deg)")

    axes[0, 1].plot(
        wheel_travel_z,
        toe_deg,
        marker="o",
        markersize=4,
    )

    axes[0, 1].scatter(
        [design_wheel_travel_z],
        [design_values["toe_deg"]],
        color="red",
        s=70,
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
    )
    annotate_design_value(
        axes[0, 1],
        design_wheel_travel_z,
        design_values["toe_deg"],
        f"{design_values['toe_deg']:.2f}°",
    )

    axes[0, 1].set_title("Toe vs Wheel Travel")

    axes[0, 1].set_ylabel("Toe(deg)")

    # --------------------------------------------------------
    # Row 2: Caster and KPI
    # --------------------------------------------------------

    axes[1, 0].plot(
        wheel_travel_z,
        caster_deg,
        marker="o",
        markersize=4,
    )

    axes[1, 0].scatter(
        [design_wheel_travel_z],
        [design_values["caster_deg"]],
        color="red",
        s=70,
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
    )
    annotate_design_value(
        axes[1, 0],
        design_wheel_travel_z,
        design_values["caster_deg"],
        f"{design_values['caster_deg']:.2f}°",
    )

    axes[1, 0].set_title("Caster vs Wheel Travel")

    axes[1, 0].set_ylabel("Caster (deg)")

    axes[1, 1].plot(
        wheel_travel_z,
        kpi_deg,
        marker="o",
        markersize=4,
    )

    axes[1, 1].scatter(
        [design_wheel_travel_z],
        [design_values["kpi_deg"]],
        color="red",
        s=70,
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
    )
    annotate_design_value(
        axes[1, 1],
        design_wheel_travel_z,
        design_values["kpi_deg"],
        f"{design_values['kpi_deg']:.2f}°",
    )

    axes[1, 1].set_title("KPI vs Wheel Travel")

    axes[1, 1].set_ylabel("KPI (deg)")

    # --------------------------------------------------------
    # Row 3: Scrub radius and mechanical trail
    # --------------------------------------------------------

    axes[2, 0].plot(
        wheel_travel_z,
        scrub_radius_mm,
        marker="o",
        markersize=4,
    )

    axes[2, 0].scatter(
        [design_wheel_travel_z],
        [design_values["scrub_radius_mm"]],
        color="red",
        s=70,
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
    )
    annotate_design_value(
        axes[2, 0],
        design_wheel_travel_z,
        design_values["scrub_radius_mm"],
        f"{design_values['scrub_radius_mm']:.2f}",
    )

    axes[2, 0].set_title(
        "Scrub Radius vs Wheel Travel"
    )

    axes[2, 0].set_xlabel(
        "Vertical wheel travel (mm)"
    )

    axes[2, 0].set_ylabel(
        "Scrub radius (mm)"
    )

    axes[2, 1].plot(
        wheel_travel_z,
        mechanical_trail_mm,
        marker="o",
        markersize=4,
    )

    axes[2, 1].scatter(
        [design_wheel_travel_z],
        [design_values["mechanical_trail_mm"]],
        color="red",
        s=70,
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
    )
    annotate_design_value(
        axes[2, 1],
        design_wheel_travel_z,
        design_values["mechanical_trail_mm"],
        f"{design_values['mechanical_trail_mm']:.2f}",
    )

    axes[2, 1].set_title(
        "Mechanical Trail vs Wheel Travel"
    )

    axes[2, 1].set_xlabel(
        "Vertical wheel travel (mm)"
    )

    axes[2, 1].set_ylabel(
        "Mechanical trail (mm)"
    )

    # --------------------------------------------------------
    # Apply common formatting to every subplot
    # --------------------------------------------------------

    for ax in axes.flat:
        ax.axvline(
            x=0.0,
            color="black",
            linestyle="--",
            linewidth=0.8,
        )

        ax.grid(True)

    fig.suptitle(
        "Double-Wishbone Suspension Kinematics",
        fontsize=16,
    )

    fig.tight_layout(
        rect=[0.0, 0.0, 1.0, 0.97]
    )

    return fig, axes

# ============================================================
# 26. PLOT THE 3D SUSPENSION ON ONE SIDE
# ============================================================

def plot_suspension_side(ax, model, results, side_name, show_labels=True):

    lower_ball_joint = results["lower_ball_joint"]
    upper_ball_joint = results["upper_ball_joint"]
    tie_rod_outer = results["tie_rod_outer"]
    tie_rod_inner = model["tie_rod_inner"]
    wheel_center = results["wheel_center"]
    hub_axis_point = results["hub_axis_point"]

    chassis_points = np.array([
        model["LCA_front_inner"],
        model["LCA_rear_inner"],
        model["UCA_front_inner"],
        model["UCA_rear_inner"],
    ])

    ball_joints = np.array([
        lower_ball_joint,
        upper_ball_joint,
    ])

    tie_rod_ends = np.array([
        tie_rod_inner,
        tie_rod_outer,
    ])

    # Only label plotted objects when requested, preventing duplicate legend entries.
    chassis_label = (
        f"{side_name} chassis pickups"
        if show_labels
        else None
    )

    ball_joint_label = (
        f"{side_name} ball joints"
        if show_labels
        else None
    )

    tie_rod_end_label = (
        f"{side_name} tie rod ends"
        if show_labels
        else None
    )

    wheel_center_label = (
        f"{side_name} wheel centre"
        if show_labels
        else None
    )

    # --------------------------------------------------------
    # Chassis pickups
    # --------------------------------------------------------

    ax.scatter(
        chassis_points[:, 0],
        chassis_points[:, 1],
        chassis_points[:, 2],
        color="blue",
        s=60,
        label=chassis_label,
    )

    # --------------------------------------------------------
    # Ball joints
    # --------------------------------------------------------

    ax.scatter(
        ball_joints[:, 0],
        ball_joints[:, 1],
        ball_joints[:, 2],
        color="red",
        s=70,
        label=ball_joint_label,
    )

    # --------------------------------------------------------
    # Tie-rod ends
    # --------------------------------------------------------

    ax.scatter(
        tie_rod_ends[:, 0],
        tie_rod_ends[:, 1],
        tie_rod_ends[:, 2],
        color="orange",
        s=70,
        label=tie_rod_end_label,
    )

    # --------------------------------------------------------
    # Wheel centre
    # --------------------------------------------------------

    ax.scatter(
        wheel_center[0],
        wheel_center[1],
        wheel_center[2],
        color="cyan",
        edgecolor="black",
        s=100,
        label=wheel_center_label,
    )

    # --------------------------------------------------------
    # Lower control arm
    # --------------------------------------------------------

    for inner_point in [
        model["LCA_front_inner"],
        model["LCA_rear_inner"],
    ]:
        ax.plot(
            [inner_point[0], lower_ball_joint[0]],
            [inner_point[1], lower_ball_joint[1]],
            [inner_point[2],lower_ball_joint[2]],
            color="green",
            linewidth=2.5,
        )

    # Lower-arm inboard pivot axis
    ax.plot(
        [model["LCA_front_inner"][0], model["LCA_rear_inner"][0]],
        [model["LCA_front_inner"][1], model["LCA_rear_inner"][1]],
        [model["LCA_front_inner"][2], model["LCA_rear_inner"][2]],
        color="green",
        linewidth=2.5,
        label=(
            f"{side_name} lower control arm"
            if show_labels
            else None
        ),
    )

    # --------------------------------------------------------
    # Upper control arm
    # --------------------------------------------------------

    for inner_point in [
        model["UCA_front_inner"],
        model["UCA_rear_inner"],
    ]:
        ax.plot(
            [inner_point[0], upper_ball_joint[0]],
            [inner_point[1], upper_ball_joint[1]],
            [inner_point[2], upper_ball_joint[2]],
            color="magenta",
            linewidth=2.5,
        )

    # Upper-arm inboard pivot axis
    ax.plot(
        [model["UCA_front_inner"][0], model["UCA_rear_inner"][0]],
        [model["UCA_front_inner"][1], model["UCA_rear_inner"][1]],
        [model["UCA_front_inner"][2], model["UCA_rear_inner"][2]],
        color="magenta",
        linewidth=2.5,
        label=(
            f"{side_name} upper control arm"
            if show_labels
            else None
        ),
    )

    # --------------------------------------------------------
    # Upright / steering axis
    # --------------------------------------------------------

    ax.plot(
        [lower_ball_joint[0], upper_ball_joint[0]],
        [lower_ball_joint[1], upper_ball_joint[1]],
        [lower_ball_joint[2], upper_ball_joint[2]],
        color="red",
        linewidth=3.0,
        label=(
            f"{side_name} steering axis"
            if show_labels
            else None
        ),
    )

    # --------------------------------------------------------
    # Tie rod
    # --------------------------------------------------------

    ax.plot(
        [tie_rod_inner[0], tie_rod_outer[0]],
        [tie_rod_inner[1], tie_rod_outer[1]],
        [tie_rod_inner[2], tie_rod_outer[2]],
        color="orange",
        linewidth=3.0,
        label=(
            f"{side_name} tie rod"
            if show_labels
            else None
        ),
    )

    # --------------------------------------------------------
    # Wheel rotation axis
    # --------------------------------------------------------

    ax.plot(
        [wheel_center[0], hub_axis_point[0]],
        [wheel_center[1], hub_axis_point[1]],
        [wheel_center[2], hub_axis_point[2]],
        color="gold",
        linewidth=3.0,
        label=(
            f"{side_name} wheel axis"
            if show_labels
            else None
        ),
    )

# ============================================================
# 27. COLLECT ALL POINTS FROM BOTH SIDES OF THE SUSPENSION
# ============================================================
def collect_axle_plot_points(
    left_model,
    left_results,
    right_model,
    right_results,
):

    model_point_names = [
        "LCA_front_inner",
        "LCA_rear_inner",
        "UCA_front_inner",
        "UCA_rear_inner",
        "tie_rod_inner",
    ]

    result_point_names = [
        "lower_ball_joint",
        "upper_ball_joint",
        "tie_rod_outer",
        "wheel_center",
        "hub_axis_point",
    ]

    points = []

    for model in [left_model, right_model]:
        for name in model_point_names:
            points.append(model[name])

    for results in [left_results, right_results]:
        for name in result_point_names:
            points.append(results[name])

    return np.array(points)

# ============================================================
# 28. ADD FUNCTION FOR EQUAL AXIS SCALING
# ============================================================

def set_3d_axes_equal(ax, points):

    coordinate_minimum = np.min(points, axis=0)

    coordinate_maximum = np.max(points, axis=0)

    coordinate_midpoint = (coordinate_minimum + coordinate_maximum) / 2.0

    coordinate_range = (coordinate_maximum - coordinate_minimum)

    maximum_range = np.max(coordinate_range)

    half_range = maximum_range / 2.0

    ax.set_xlim(
        coordinate_midpoint[0] - half_range,
        coordinate_midpoint[0] + half_range,
    )

    ax.set_ylim(
        coordinate_midpoint[1] - half_range,
        coordinate_midpoint[1] + half_range,
    )

    ax.set_zlim(
        coordinate_midpoint[2] - half_range,
        coordinate_midpoint[2] + half_range,
    )

    try:
        ax.set_box_aspect(
            [1.0, 1.0, 1.0]
        )
    except AttributeError:
        pass

# ============================================================
# 29. PLOT WHEELS
# ============================================================

def plot_wheel_circle(ax, wheel_center, hub_axis_point, wheel_radius, color="black"):

    wheel_axis = unit_vector(hub_axis_point - wheel_center)

    global_up = np.array([0.0, 0.0, 1.0])

    circle_axis_1 = global_up - (np.dot(global_up, wheel_axis)* wheel_axis)

    if vector_length(circle_axis_1) < 1e-12:
        alternative = np.array([
            1.0,
            0.0,
            0.0,
        ])

        circle_axis_1 = alternative - (np.dot(alternative, wheel_axis)* wheel_axis)

    circle_axis_1 = unit_vector(circle_axis_1)

    circle_axis_2 = unit_vector(np.cross(wheel_axis, circle_axis_1))

    theta = np.linspace(0.0, 2.0 * np.pi, 100)

    wheel_points = []

    for angle in theta:
        point = (
            wheel_center
            + wheel_radius
            * np.cos(angle)
            * circle_axis_1
            + wheel_radius
            * np.sin(angle)
            * circle_axis_2
        )

        wheel_points.append(point)

    wheel_points = np.array(wheel_points)

    ax.plot(
        wheel_points[:, 0],
        wheel_points[:, 1],
        wheel_points[:, 2],
        color=color,
        linewidth=2.0,
    )

# ============================================================
# 30. PLOT THE 3D SUSPENSION ON BOTH SIDES
# ============================================================

def plot_suspension_axle_3d(left_model, left_results, right_model, right_results):

    fig = plt.figure(
        figsize=(13, 9)
    )

    ax = fig.add_subplot(
        111,
        projection="3d",
    )

    # Plot left suspension
    plot_suspension_side(ax, left_model, left_results, "Left", True)

    # Plot right suspension
    plot_suspension_side(ax, right_model, right_results, "Right", False)

    # Plot left wheel
    plot_wheel_circle(ax, left_results["wheel_center"], left_results["hub_axis_point"], left_model["wheel_radius_mm"])

    #Plot right wheel
    plot_wheel_circle(ax, right_results["wheel_center"], right_results["hub_axis_point"], right_model["wheel_radius_mm"])

    # --------------------------------------------------------
    # Connect corresponding chassis pickups
    # --------------------------------------------------------

    chassis_point_names = [
        "LCA_front_inner",
        "LCA_rear_inner",
        "UCA_front_inner",
        "UCA_rear_inner",
        "tie_rod_inner",
    ]

    for point_name in chassis_point_names:
        left_point = left_model[point_name]
        right_point = right_model[point_name]

        ax.plot(
            [left_point[0], right_point[0]],
            [left_point[1], right_point[1]],
            [left_point[2], right_point[2]],
            color="gray",
            linestyle="--",
            linewidth=1.0,
            alpha=0.6,
        )

    # --------------------------------------------------------
    # Vehicle centreline
    # --------------------------------------------------------

    all_points = collect_axle_plot_points(
        left_model,
        left_results,
        right_model,
        right_results,
    )

    minimum_z = np.min(all_points[:, 2])
    maximum_z = np.max(all_points[:, 2])

    ax.plot(
        [0.0, 0.0],
        [0.0, 0.0],
        [minimum_z, maximum_z],
        color="black",
        linestyle="--",
        linewidth=1.5,
        label="Vehicle centre plane",
    )

    # --------------------------------------------------------
    # Labels and title
    # --------------------------------------------------------

    ax.set_xlabel("X, forward (mm)")
    ax.set_ylabel("Y, left (mm)")
    ax.set_zlabel("Z, upward (mm)")

    ax.set_title(
        "Complete Double-Wishbone Front Axle\n"
        f"Left LCA angle = "
        f"{left_results['lca_angle_deg']:.1f} deg, "
        f"Right LCA angle = "
        f"{right_results['lca_angle_deg']:.1f} deg"
    )

    # --------------------------------------------------------
    # Set equal axis scaling
    # --------------------------------------------------------

    set_3d_axes_equal(
        ax,
        all_points,
    )

    ax.legend(
        loc="upper left",
        fontsize=8,
    )

    ax.grid(True)

    fig.tight_layout()

    return fig, ax


# ============================================================
# 31. PLOT FRONT VIEW ROLL CENTER
# ============================================================
def plot_front_view_roll_center(left_model, left_results, right_model, right_results, roll_center_yz, construction):

    # Read suspension points and convert them to front view
    left_lbj_yz = point_to_yz(left_results["lower_ball_joint"])
    left_ubj_yz = point_to_yz(left_results["upper_ball_joint"])
    right_lbj_yz = point_to_yz(right_results["lower_ball_joint"])
    right_ubj_yz = point_to_yz(right_results["upper_ball_joint"])

    left_lca_rear_yz = point_to_yz(left_model["LCA_rear_inner"])
    left_lca_front_yz = point_to_yz(left_model["LCA_front_inner"])
    right_lca_rear_yz = point_to_yz(right_model["LCA_rear_inner"])
    right_lca_front_yz = point_to_yz(right_model["LCA_front_inner"])

    left_uca_rear_yz = point_to_yz(left_model["UCA_rear_inner"])
    left_uca_front_yz = point_to_yz(left_model["UCA_front_inner"])
    right_uca_rear_yz = point_to_yz(right_model["UCA_rear_inner"])
    right_uca_front_yz = point_to_yz(right_model["UCA_front_inner"])

    # Read roll-center construction points
    left_ic_yz = construction["left_instant_center_yz"]
    right_ic_yz = construction["right_instant_center_yz"]
    left_contact_yz = construction["left_contact_patch_yz"]
    right_contact_yz = construction["right_contact_patch_yz"]

    # Determine ground coordinate
    common_ground_z = 0.5 * (left_contact_yz[1]+ right_contact_yz[1])
    roll_center_height_mm = roll_center_yz[1]- common_ground_z

    # Create the figure
    fig, ax = plt.subplots(figsize=(13, 8))

    # Plot ground and vehicle centerline
    ax.axhline(
        y=common_ground_z,
        color="black",
        linestyle="-",
        linewidth=1.5,
        label="Ground",
    )

    ax.axvline(
        x=0.0,
        color="gray",
        linestyle="--",
        linewidth=1.2,
        label="Vehicle centerline",
    )

    # Plot ball joints
    ax.scatter(
        [left_lbj_yz[0], right_lbj_yz[0]],
        [left_lbj_yz[1], right_lbj_yz[1]],
        color="green",
        edgecolor="black",
        s=75,
        zorder=5,
        label="Lower ball joints",
    )

    ax.scatter(
        [left_ubj_yz[0], right_ubj_yz[0]],
        [left_ubj_yz[1], right_ubj_yz[1]],
        color="magenta",
        edgecolor="black",
        s=75,
        zorder=5,
        label="Upper ball joints",
    )

    # Plot front inner pick up points
    ax.scatter(
        [left_lca_front_yz[0], right_lca_front_yz[0], left_uca_front_yz[0], right_uca_front_yz[0]],
        [left_lca_front_yz[1], right_lca_front_yz[1], left_uca_front_yz[1], right_uca_front_yz[1]],
        color="blue",
        edgecolor="black",
        s=75,
        zorder=3,
        label="Inner pick up points",
    )

    ax.scatter(
        [left_lca_rear_yz[0], right_lca_rear_yz[0], left_uca_rear_yz[0], right_uca_rear_yz[0]],
        [left_lca_rear_yz[1], right_lca_rear_yz[1], left_uca_rear_yz[1], right_uca_rear_yz[1]],
        color="blue",
        edgecolor="black",
        s=75,
        zorder=2,
        alpha = 0.2
    )

    # Plot the left and right uprights
    ax.plot(
        [left_lbj_yz[0], left_ubj_yz[0]],
        [left_lbj_yz[1], left_ubj_yz[1]],
        color="red",
        linewidth=2.5,
        label="Uprights",
    )

    ax.plot(
        [right_lbj_yz[0], right_ubj_yz[0]],
        [right_lbj_yz[1], right_ubj_yz[1]],
        color="red",
        linewidth=2.5,
    )

    # Plot contact points
    ax.scatter(
        [left_contact_yz[0], right_contact_yz[0]],
        [left_contact_yz[1], right_contact_yz[1]],
        color="black",
        marker="s",
        s=75,
        zorder=6,
        label="Contact points",
    )

    # Plot instant centers
    ax.scatter(
        [left_ic_yz[0], right_ic_yz[0]],
        [left_ic_yz[1], right_ic_yz[1]],
        color="orange",
        edgecolor="black",
        marker="X",
        s=130,
        zorder=7,
        label="Instant centers",
    )

    # Plot force lines
    ax.plot(
        [left_contact_yz[0], left_ic_yz[0]],
        [left_contact_yz[1], left_ic_yz[1]],
        color="dodgerblue",
        linestyle="--",
        linewidth=2.0,
        label="Contact point to instant-center lines",
    )

    ax.plot(
        [right_contact_yz[0], right_ic_yz[0]],
        [right_contact_yz[1], right_ic_yz[1]],
        color="dodgerblue",
        linestyle="--",
        linewidth=2.0,
    )

    # Plot ball joints to IC lines
    ax.plot(
        [left_lbj_yz[0], left_ic_yz[0]],
        [left_lbj_yz[1], left_ic_yz[1]],
        color="green",
        linestyle="--",
        linewidth=2.0,
        label="Ball joints to instant-center lines",
    )

    ax.plot(
        [left_ubj_yz[0], left_ic_yz[0]],
        [left_ubj_yz[1], left_ic_yz[1]],
        color="green",
        linestyle="--",
        linewidth=2.0,
    )

    ax.plot(
        [right_lbj_yz[0], right_ic_yz[0]],
        [right_lbj_yz[1], right_ic_yz[1]],
        color="green",
        linestyle="--",
        linewidth=2.0, 
    )

    ax.plot(
        [right_ubj_yz[0], right_ic_yz[0]],
        [right_ubj_yz[1], right_ic_yz[1]],
        color="green",
        linestyle="--",
        linewidth=2.0, 
    )

    # Plot roll center
    ax.scatter(
        roll_center_yz[0],
        roll_center_yz[1],
        color="red",
        edgecolor="black",
        marker="*",
        s=260,
        zorder=10,
        label="Front roll center",
    )

    # Annotate important points
    ax.annotate(
        "Left IC",
        xy=left_ic_yz,
        xytext=(-40, 8),
        textcoords="offset points",
    )

    ax.annotate(
        "Right IC",
        xy=right_ic_yz,
        xytext=(8, 8),
        textcoords="offset points",
    )

    ax.annotate(
        "Roll center\n"
        f"Y = {roll_center_yz[0]:.2f} mm\n"
        f"Z = {roll_center_yz[1]:.2f} mm\n"
        f"Height = {roll_center_height_mm:.2f} mm",
        xy=roll_center_yz,
        xytext=(-50, -70),
        textcoords="offset points",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "edgecolor": "red",
            "alpha": 0.9,
        },
        arrowprops={
            "arrowstyle": "->",
            "color": "red",
        },
    )

    # Plot formatting
    ax.set_xlabel("Global Y coordinate, left positive (mm)")
    ax.set_ylabel("Global Z coordinate, upward positive (mm)")
    ax.set_title(
        "Front-View Roll-Center Construction\n"
        f"Roll-center height above ground = "
        f"{roll_center_height_mm:.2f} mm"
    )

    ax.grid(True)

    # Equal Y-Z scaling so geometric angles are not distorted.
    ax.set_aspect(
        "equal",
        adjustable="datalim",
    )

    ax.legend(
        loc="best",
        fontsize=9,
    )

    fig.tight_layout()

    return fig, ax

# ============================================================
# 32. PLOT ROLL CENTER GRAPHS
# ============================================================

def plot_roll_center_height(roll_center_sweep):

    wheel_travel_z = np.array([
        result["average_wheel_travel_z_mm"]
        for result in roll_center_sweep
    ])

    roll_center_height = np.array([
        result["roll_center_height_mm"]
        for result in roll_center_sweep
    ])

    # Sort data by wheel travel
    sort_order = np.argsort(wheel_travel_z)

    wheel_travel_z = wheel_travel_z[sort_order]
    roll_center_height = roll_center_height[sort_order]

    # Find design position
    design_index = np.argmin(
        np.abs(wheel_travel_z)
    )

    design_wheel_travel = wheel_travel_z[design_index]
    design_roll_center_height = roll_center_height[design_index]

    # Create figure
    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    # Plot roll-center migration
    ax.plot(
        wheel_travel_z,
        roll_center_height,
        marker="o",
        markersize=4,
    )

    # Highlight design position
    ax.scatter(
        design_wheel_travel,
        design_roll_center_height,
        color="red",
        edgecolor="black",
        s=90,
        zorder=10,
    )

    # Annotate design roll-center height
    ax.annotate(
        f"{design_roll_center_height:.1f} mm",
        xy=(
            design_wheel_travel,
            design_roll_center_height,
        ),
        xytext=(10, 10),
        textcoords="offset points",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "edgecolor": "red",
            "alpha": 0.9,
        },
    )

    # Mark design wheel position
    ax.axvline(
        x=0.0,
        color="black",
        linestyle="--",
        linewidth=0.8,
    )

    ax.set_title(
        "Front Roll-Center Height vs Wheel Travel"
    )

    ax.set_xlabel(
        "Vertical wheel travel (mm)"
    )

    ax.set_ylabel(
        "Roll-center height above ground (mm)"
    )

    ax.grid(True)

    fig.tight_layout()

    return fig, ax

# ============================================================
# 33. PRINT MODEL INFORMATION
# ============================================================

def print_model_information(model):
    """
    Print the rigid-link lengths calculated from the input points.
    """

    print("\nCalculated design link lengths")
    print("--------------------------------")

    print(
        f"LCA front leg: "
        f"{model['LCA_front_length']:.3f} mm"
    )

    print(
        f"LCA rear leg:  "
        f"{model['LCA_rear_length']:.3f} mm"
    )

    print(
        f"UCA front leg: "
        f"{model['UCA_front_length']:.3f} mm"
    )

    print(
        f"UCA rear leg:  "
        f"{model['UCA_rear_length']:.3f} mm"
    )

    print(
        f"Tie rod:  "
        f"{model['tie_rod_length']:.3f} mm"
    )

    print(
        f"Upright length: "
        f"{model['upright_length']:.3f} mm"
    )


def print_constraint_check(model, results):
    """
    Print the maximum constraint error for one solved position.
    """

    errors = calculate_constraint_errors(model, results)

    print(
        f"\nConstraint check at "
        f"{results['lca_angle_deg']:.2f} deg"
    )
    print("--------------------------------")

    for error_name, error_value in errors.items():
        print(
            f"{error_name:22s}: "
            f"{error_value:+.9f} mm"
        )

    maximum_error = max(
        abs(error_value)
        for error_value in errors.values()
    )

    print(
        f"Maximum absolute error: "
        f"{maximum_error:.9f} mm"
    )


# ============================================================
# 34. MAIN PROGRAM
# ============================================================

def main():
    """
    Main program execution.
    """
    os.makedirs(
        OUTPUT_FOLDER,
        exist_ok=True,
    )

    left_model = create_suspension_model(hardpoints)

    right_hardpoints = mirror_hardpoints(hardpoints)
    right_model = create_suspension_model(right_hardpoints)

    print_model_information(left_model)

    # Generate the requested LCA angle values
    lca_angles_deg = np.linspace(
        ANGLE_START_DEG,
        ANGLE_END_DEG,
        NUMBER_OF_POSITIONS,
    )

    # Generate angles for both sides of suspension
    left_angles_deg = lca_angles_deg
    right_angles_deg = -lca_angles_deg

    # Generate DISPLAY angles for both sides of suspension
    left_display_angle_deg = DISPLAY_ANGLE_DEG
    right_display_angle_deg = -DISPLAY_ANGLE_DEG

    # Solve all suspension positions
    left_sweep_results = run_angle_sweep(left_model,left_angles_deg)

    right_sweep_results = run_angle_sweep(right_model, right_angles_deg)
    right_sweep_results.reverse()

    # Calculate roll center for sweeped positions
    roll_center_sweep = calculate_roll_center_sweep(left_sweep_results,right_sweep_results,left_model,right_model)

    # Solve the positions used for the 3D display
    left_display_results = solve_suspension_position(
        left_model,
        left_display_angle_deg,
        left_model["upper_ball_joint"],
        left_model["tie_rod_outer"],
    )
    left_display_results = add_derived_geometry(left_model, left_display_results)

    right_display_results = solve_suspension_position(
        right_model,
        right_display_angle_deg,
        right_model["upper_ball_joint"],
        right_model["tie_rod_outer"],
    )
    right_display_results = add_derived_geometry(right_model, right_display_results)

    # Calculate roll center for display
    front_roll_center_yz, roll_center_construction = (
        calculate_roll_center(
            left_model,
            left_display_results,
            right_model,
            right_display_results,
        )
    )

    # Check that rigid-link constraints are maintained
    print_constraint_check(
        left_model,
        left_display_results,
    )

    # Produce plots
    # --------------------------------------------------------
    # 1. Suspension kinematic curves
    # --------------------------------------------------------

    fig_kinematics, axes_kinematics = plot_kinematic_curves(
        right_model,
        right_sweep_results,
    )

    fig_kinematics.savefig(
        os.path.join(
            OUTPUT_FOLDER,
            "suspension_kinematic_curves.png",
        ),
        dpi=FIGURE_DPI,
        bbox_inches="tight",
    )

    # --------------------------------------------------------
    # 2. 3D suspension model
    # --------------------------------------------------------

    fig_3d, ax_3d = plot_suspension_axle_3d(
        left_model,
        left_display_results,
        right_model,
        right_display_results,
    )

    fig_3d.savefig(
        os.path.join(
            OUTPUT_FOLDER,
            "suspension_3d_model.png",
        ),
        dpi=FIGURE_DPI,
        bbox_inches="tight",
    )

    # --------------------------------------------------------
    # 3. Front-view roll-centre construction
    # --------------------------------------------------------

    fig_roll_center, ax_roll_center = plot_front_view_roll_center(
        left_model,
        left_display_results,
        right_model,
        right_display_results,
        front_roll_center_yz,
        roll_center_construction,
    )

    fig_roll_center.savefig(
        os.path.join(
            OUTPUT_FOLDER,
            "front_roll_center_construction.png",
        ),
        dpi=FIGURE_DPI,
        bbox_inches="tight",
    )

    # --------------------------------------------------------
    # 4. Roll-centre migration
    # --------------------------------------------------------

    fig_rc_migration, axes_rc_migration = plot_roll_center_height(roll_center_sweep)

    fig_rc_migration.savefig(
        os.path.join(
            OUTPUT_FOLDER,
            "roll_center_migration.png",
        ),
        dpi=FIGURE_DPI,
        bbox_inches="tight",
    )


    # Confirm that figures have been saved
    print("\nSaved portfolio figures:")
    print("  suspension_kinematic_curves.png")
    print("  suspension_3d_model.png")
    print("  front_roll_center_construction.png")
    print("  roll_center_migration.png")


    # Display all figures after saving
    plt.show()


    # Check mirrored positions
    check_mirrored_positions(left_display_results,right_display_results,)

if __name__ == "__main__":
    main()