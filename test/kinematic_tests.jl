@testset "KinematicTests" begin
    # testing 100 random poses in direct-inverse-direct fashion on outer kinematics
    old_man = build_old_manipulator()
    for i = 1:100
        C = forward_kinematics(old_man, π .- 2π*rand(5))
        θ1 = inverse_kinematics(old_man, C, sol = 1)
        θ2 = inverse_kinematics(old_man, C, sol = 2)
        @test isapprox(C.mat, forward_kinematics(old_man, θ1).mat, atol = 1e-6) || isapprox(C.mat, forward_kinematics(old_man, θ2).mat, atol = 1e-6)
    end
end