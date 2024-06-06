using InstantFrame

# Define the material properties:
Material = InstantFrame.Material(
    names = ["Steel"], 
    E     = [29000], 
    ν     = [0.3], 
    ρ     = [0])

# Define the cross-sectional properties:
CrossSection = InstantFrame.CrossSection(
    names = ["Beam"], 
    A     = [9.12], 
    Iy    = [0], 
    Iz    = [110], 
    J     = [0])

# Define the element-to-element connections:
Connection = InstantFrame.Connection(
    names     = ["Rigid"], 
    stiffness = (ux = [Inf], uy = [Inf], uz = [Inf], rx = [Inf], ry = [Inf], rz = [Inf]))

# Define the nodes:
Node = InstantFrame.Node(
    numbers     = 1:11, 
    coordinates = [
        ( 0 * 18, 0, 0), 
        ( 1 * 18, 0, 0), 
        ( 2 * 18, 0, 0), 
        ( 3 * 18, 0, 0), 
        ( 4 * 18, 0, 0), 
        ( 5 * 18, 0, 0), 
        ( 6 * 18, 0, 0), 
        ( 7 * 18, 0, 0), 
        ( 8 * 18, 0, 0), 
        ( 9 * 18, 0, 0), 
        (10 * 18, 0, 0)])

# Define the elements:
Element = InstantFrame.Element(
    types         = ["", "", "", "", "", "", "", "", "", ""],
    numbers       = 1:10, 
    nodes         = [
        ( 1,  2), 
        ( 2,  3), 
        ( 3,  4), 
        ( 4,  5), 
        ( 5,  6), 
        ( 6,  7), 
        ( 7,  8), 
        ( 8,  9), 
        ( 9, 10), 
        (10, 11)], 
    orientation   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    connections   = [
        ("Rigid", "Rigid"), 
        ("Rigid", "Rigid"), 
        ("Rigid", "Rigid"), 
        ("Rigid", "Rigid"), 
        ("Rigid", "Rigid"), 
        ("Rigid", "Rigid"), 
        ("Rigid", "Rigid"), 
        ("Rigid", "Rigid"), 
        ("Rigid", "Rigid"), 
        ("Rigid", "Rigid")], 
    cross_section = [ "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam"], 
    material      = ["Steel", "Steel", "Steel", "Steel", "Steel", "Steel", "Steel", "Steel", "Steel", "Steel"])

# Define the boundary conditions:
Support = InstantFrame.Support(
    nodes     = 1:11, 
    stiffness = (
        uX = [Inf,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0], 
        uY = [Inf,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0], 
        uZ = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf],
        rX = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf], 
        rY = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf], 
        rZ = [Inf,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0]))

# Define the distributed loads:
UniformLoad = InstantFrame.UniformLoad(nothing)

# Define the concentrated loads:
PointLoad = InstantFrame.PointLoad(
    labels     = [""], 
    nodes      = [11], 
    magnitudes = (FX = [-50], FY = [-1], FZ = [0], MX = [0], MY = [0], MZ = [0]))

# Solve:
Model = InstantFrame.solve(Node, CrossSection, Material, Connection, Element, Support, UniformLoad, PointLoad, analysis_type = "second order", solution_tolerance = 1E-6)

# Get the vertical displacement at the free end:
Δ = -Model.solution.displacements[11][2]