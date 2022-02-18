module InstantFrame

struct Node

    x::Float64
    y::Float64
    z::Float64

end

struct Beam

    start_node::Int64
    end_node::Int64
    start_connection::Connection
    end_connection::Connection
    cross_section::String
    material::String

end

struct Connection

    rx::Float64
    ry::Float64
    rz::Float64 

end

struct Material

    name::String
    E::Float64
    ν::Float64

end

struct CrossSection

    name::String
    A::Float64
    Ixx::Float64
    Iyy::Float64
    J::Float64

end

#calculate element orientations

#define element elastic stiffness matrices

#define element geometric stiffness matrices

#define element mass matrices

#apply member end conditions

#calculate element global stiffness matrices

#assemble Ke and Kg

#assemble M

#apply boundary conditions

#solve for displacement field

#solve for internal forces

#solve for vibration modes

end # module
