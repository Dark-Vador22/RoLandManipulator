# altering the relative paths in URDF-file, once this package is loaded
filename = abspath(joinpath(@__DIR__, "..", "urdf/SHI-100-000-000_ccd.urdf"))
xdoc = parse_file(filename)
xroot = LightXML.root(xdoc)
# changing attributes only of the visual - there is also collision!
for link in get_elements_by_tagname(xroot, "link")
    mesh = get_elements_by_tagname(get_elements_by_tagname(get_elements_by_tagname(link, "visual")[1], "geometry")[1], "mesh")[1]
    value = attribute(mesh, "filename")[2:end]
    ### NEVER USE JOINPATH HERE!!!!
    new_value = abspath(String(@__DIR__)*"/../"*value)
    set_attribute(mesh, "filename", new_value)
end
# save file as copy in package
if Sys.islinux() || Sys.isapple()
    save_file(xdoc, joinpath(@__DIR__, "..", "urdf/local_copy.urdf"))
elseif Sys.iswindows()
    save_file(xdoc, joinpath(@__DIR__, "..", "urdf\\local_copy.urdf"))
else
    @warn "Can't write local URDF copy - visualisation will not work"
end

"""
Definition of joints and links by Lie group representation
"""
function outer_geometry()
     # defining geometry and putting URDF path
     C0 = SE3(zeros(6))

     # defining joints by screws
     joint_axis = [se3([0, 0, 0, 1, 0, 0]),
     se3([1, 0, 0, 0, 0, 0]),
     se3([1, 0, 0, 0, 0, 0]),
     se3([1, 0, 0, 0, 0, 0]),
     se3([1, 0, 0, 0, 0, 0])]
 
     # defining relative transforms
     Cl = [SE3([0, 0, 0, -0.5, 0.5, 0.05]),
     SE3(zeros(6)),
     SE3([0, 0, 0, 0, 0, -0.4]),
     SE3([0, 0, 0, 0, 0, -0.4]),
     SE3([0, π/2, 0, 0, 0, 0]),
     SE3([0, 0, 0, 0.1, 0, 0])]
     
     return C0, Cl, joint_axis
end

"""
Builds the old manipulator (inner belt structure)
"""
function build_manipulator()
    
    # general joint space quantities
    C0, Cl, joint_axis = outer_geometry()   

    # teeth numbers of gear stages
    n1 = 16
    n2 = 16
    n3 = 24
    n4 = 24
    n5i = 48
    n5o = 24
    n6 = 48
    n7i = 48
    n7o = 28
    n8 = 42
    n9i = 24
    n9o = 20
    n10i = 24
    n10o = 20
    n11 = 20
    n12 = 20
    
    # building structure matrix
    g1 = (n1/n5i)*(n5o/n6)
    g2 = (n2/n7i)*(n7o/n8)
    g3 = n3/n9i
    g4 = n4/n10i
    g5 = (n3/n9i)*(n9o/n11)
    g6 = (n4/n10i)*(n10o/n12)
    
    S = [[1/g1, 1, 1, 1] [0, 1/g2, 1/g3, 1/g4] [0, 0, 1/g5, 1/g6] [0, 0, 1/g5, -1/g6]]
    
    # building actuation space stiffness matrix for planar 4 DOFs
    # length specific stiffness of 10mm and 16mm belt
    ks_10 = 787/4e-3
    ks_16 = 1342/4e-3
    
    # length in belt stages
    l1_5 = 57e-3                    # 10
    l5_6 = 68e-3                    # 16
    l2_7 = 57e-3                    # 10
    l7_8 = (315e-3 + 323e-3)/2      # 16
    l3_9 = (398e-3 + 400e-3)/2      # 10
    l4_10 = (398e-3 + 400e-3)/2     # 10
    l9_11 = (395e-3 + 398e-3)/2     # 10
    l10_12 = (395e-3 + 398e-3)/2    # 10
    
    # radius from teeth number with module 5mm
    rs = [n1, n2, n3, n4]*5e-3/2π

    # belt stiffness in all strings
    kss = [[ks_10/l1_5, ks_16/l5_6], [ks_10/l2_7, ks_16/l7_8], [ks_10/l3_9, ks_10/l9_11], [ks_10/l4_10, ks_10/l10_12]]
    
    # reduction ratios
    ρs = [[n5o/n5i], [n7o/n7i], [n9o/n9i], [n10o/n10i]]

    K_μ = diagm([actuation_space_stiffness(rs[i], kss[i], ρs[i]) for i = 1:4])

    return Manipulator(C0, Cl, joint_axis, S, 6., K_μ, abspath(joinpath(@__DIR__, "..", "urdf/local_copy.urdf")), friction = 0.1)
end
