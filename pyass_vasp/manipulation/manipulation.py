import numpy as np

def rotate_ions(lines, center, direction, angle, POSCAR='POSCAR', output='POSCAR-rotated', inplace=False):
    """
    Rotate the ions in certain lines of POSCAR around a certain point,
        along a certain direction, to a certain angle.

    Parameters
    ----------
    lines: list
        the line numbers of ions from POSCAR that you wish to rotate, counted from 1
    center: list
        the point that the rotation axis will pass, in "direct" coordinates
    direction: list
        the rotation axis' "cartesian" coordinates. No need to be normalized.
    angle: float
        the rotation angle in degree, counter-clockwise
    POSCAR: string
        POSCAR file name, default to 'POSCAR'
    output: string
        the output file name, default to 'POSCAR-rotated'
    inplace: bool
        directly change POSCAR in place without creating a new file

    Returns
    -------
    a dict, containing ion_positions_before, ion_positions_after, rotation matrix

    If you wish to use it, ion_positions_after = np.dot(ion_positions_before, A.T)
    """

    if inplace:
        output = POSCAR

    # define the rotation matrix A. Here we use the active rotation scheme, i.e. rotating the object
    [u_x, u_y, u_z] = direction / np.linalg.norm(direction)
    angle = angle / 180. * np.pi
    A = np.eye(3) * np.cos(angle) + np.array([[0, -u_z, u_y], [u_z, 0, -u_x], [-u_y, u_x, 0]]) \
                * np.sin(angle) + np.array([[u_x ** 2, u_x * u_y, u_x * u_z],
                                            [u_x * u_y, u_y ** 2, u_y * u_z],
                                            [u_x * u_z, u_y * u_z, u_z ** 2]]) * (1 - np.cos(angle))

    # np.set_printoptions(suppress=True)
    line_number_list = np.array(lines) - 1  # machine counts from 0


    # open the POSCAR that has ion positions part
    with open(POSCAR, 'r') as f:
        POSCAR = f.readlines()

    # obtain the basis vectors from POSCAR
    basis_vectors = np.zeros((3, 3))
    for i, line in enumerate(POSCAR[2:5]):
        basis_vectors[i] = line.split()

    # convert the string form of ion positions into a numpy numeral array
    ion_positions = np.zeros((len(line_number_list), 3))
    for i, line in enumerate(line_number_list):
        ion_positions[i] = POSCAR[line].split()[:3]
    # subtract the translational vector
    ion_positions -= center
    # transform the coordinates from direct to cartesian
    ion_positions = np.dot(ion_positions, basis_vectors)
    # the matrix multiplication
    ion_positions_after = np.dot(ion_positions, A.T)
    # transform the coordinates from cartesian back to direct
    ion_positions_after = np.dot(ion_positions_after, np.linalg.inv(basis_vectors))

    # ion_positions_after = ion_positions_after / np.cos(angle)

    # add the translational vector
    ion_positions_after = ion_positions_after + center

    for i, line in enumerate(line_number_list):
        POSCAR[line] = '{0[0]:11.8f}   {0[1]:11.8f}   {0[2]:11.8f}\n'.format(ion_positions_after[i])

    with open(output, 'w') as f:
        f.writelines(POSCAR)

    return {'ion_positions': ion_positions, 'ion_positions_after': ion_positions_after, 'rotation_matrix': A}


def apply_strain(lines, strain_matrix, POSCAR='POSCAR', output='POSCAR-strained', inplace=False):
    """
    Apply strain to the cell basis vectors in POSCAR. Following
        basis_vectors_after = np.dot(basis_vectors_before, A.transpose())

    Parameters
    ----------
    lines: list
        the line numbers of ions from POSCAR that you wish to rotate, counted from 1
    strain_matrix: 2D list/array
        the 2D
    POSCAR: string
        POSCAR file name, default to 'POSCAR'
    output: string
        the output file name, default to 'POSCAR-rotated'
    inplace: bool
        directly change POSCAR in place without creating a new file

    Returns
    -------
    a dict, containing basis_vectors_before, basis_vectors_after
    """

    if inplace:
        output = POSCAR

    # open the POSCAR file
    with open(POSCAR, 'r') as f:
        POSCAR = f.readlines()

    # obtain the basis vector 2D array
    basis_vectors = np.zeros((3, 3))
    for i, line in enumerate(POSCAR[2:5]):
        basis_vectors[i] = line.split()

    # apply the specified strain
    basis_vectors_after = np.dot(basis_vectors, strain_matrix)

    # replace the POSCAR variable content
    for i, line_num in enumerate(range(2, 5)):
        POSCAR[line_num] = '{0[0]:11.8f}   {0[1]:11.8f}   {0[2]:11.8f}\n'.format(basis_vectors_after[i])

    # write to it
    with open(output, 'w') as f:
        f.writelines(POSCAR)

    return {'basis_vectors_before': basis_vectors, 'basis_vectors_after': basis_vectors_after}