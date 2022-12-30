module Export 

using MAT

function mastan(node, element, filepath, filename)

    #     filepath = "/Volumes/GoogleDrive/Shared drives/RunToSolve/Projects/OneRack/production_model_validation_report/check_3_CFS_1st_order/LC1_gravity_loads/Mastan_full_model"
    #     filename = "mastan_try3.mat"
    # mastan_model = matread(joinpath(filepath, filename))
        #seed model
        mastan_model = Dict{String, Any}("sect_info" => Matrix{Float64}(undef, 0, 0), "view_settings" => [0.1 0.1 10.0 2.294558781116753 14.0 14.0 1.0 1382.0 765.0], "normincr_settings" => Matrix{Float64}(undef, 0, 0), "sec_inel_settings" => Matrix{Float64}(undef, 0, 0), "settle_info" => [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN], "bimom_settings" => Matrix{Float64}(undef, 0, 0), "mat_name" => Matrix{Float64}(undef, 0, 0), "first_el_settings" => Matrix{Float64}(undef, 0, 0), "fixity_info" => [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN], "axial_settings" => Matrix{Float64}(undef, 0, 0), "sheary_settings" => Matrix{Float64}(undef, 0, 0), "momx_settings" => Matrix{Float64}(undef, 0, 0), "elem_info" => Matrix{Float64}(undef, 0, 0), "shearz_settings" => Matrix{Float64}(undef, 0, 0), "ground_motion_data" => Matrix{Float64}(undef, 0, 0), "periods_info" => Matrix{Float64}(undef, 0, 0), "first_inel_settings" => Matrix{Float64}(undef, 0, 0), "node_info" => [0.0 0.0 0.0 113.000244140625 114.000244140625; 10.0 0.0 0.0 115.000244140625 116.000244140625; 10.0 30.0 0.0 117.000244140625 118.000244140625; 10.0 30.0 40.0 119.000244140625 120.000244140625], "deflect_settings" => Matrix{Float64}(undef, 0, 0), "ibuck_settings" => Matrix{Float64}(undef, 0, 0), "analysis_info" => Matrix{Float64}(undef, 0, 0), "nload_info" => [0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN; 0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN; 0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN; 0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN], "ebuck_settings" => Matrix{Float64}(undef, 0, 0), "momy_settings" => Matrix{Float64}(undef, 0, 0), "yldsurfval_settings" => Matrix{Float64}(undef, 0, 0), "period_settings" => Matrix{Float64}(undef, 0, 0), "mat_info" => Matrix{Float64}(undef, 0, 0), "uniload_info" => Matrix{Float64}(undef, 0, 0), "model_title" => Matrix{Float64}(undef, 0, 0), "sect_name" => Matrix{Float64}(undef, 0, 0), "first_elastic_dynamic_settings" => Matrix{Float64}(undef, 0, 0), "sec_elastic_dynamic_settings" => Matrix{Float64}(undef, 0, 0), "momz_settings" => Matrix{Float64}(undef, 0, 0), "sec_el_settings" => Matrix{Float64}(undef, 0, 0))
    
    
        num_nodes = length(node.numbers)
    
        node_matrix = Matrix{Float64}(undef, size(node.numbers, 1), 5)
    
        node_matrix[:, 1] = [node.coordinates[i][1] for i in eachindex(node.coordinates)]
        node_matrix[:, 2] = [node.coordinates[i][2] for i in eachindex(node.coordinates)]
        node_matrix[:, 3] = [node.coordinates[i][3] for i in eachindex(node.coordinates)]
        node_matrix[:, 4] = range(113.000244140625, 113.000244140625+(num_nodes-1)*2, step=2.0)
        node_matrix[:, 5] = range(114.000244140625, 114.000244140625+(num_nodes-1)*2, step=2.0)
    
        mastan_model["node_info"] = node_matrix
    
        mastan_model["fixity_info"] = fill(NaN, num_nodes, 12)
    
        mastan_model["settle_info"] = fill(NaN, num_nodes, 12)
    
        mastan_model["nload_info"] = zeros(Float64, num_nodes, 14)
        mastan_model["nload_info"][:,7:12] .= NaN 
        mastan_model["nload_info"][:,14] .= NaN
    
        # #elem_info
        # 1	2	0	0	0	561.000122070313	562.000122070313	0	1	0	563.000122070313	0	0	566.000122070313	568.000122070313	0	0	567.000122070313	569.000122070313	Inf	Inf	Inf	Inf	Inf	Inf	Inf	Inf	564.000122070313	565.000122070313
    
        num_elem = length(element.numbers)
        elem_matrix = zeros(Float64, size(element.numbers, 1), 29)
    
    
        elem_matrix[:, 1] = [element.nodes[i][1] for i in eachindex(element.nodes)]
        elem_matrix[:, 2] = [element.nodes[i][2] for i in eachindex(element.nodes)]
    
        start_index = maximum(node_matrix[:,5]) + 1
    
        # num_elem=2
        # range(start_index, start_index +(num_elem-1)*9, step=9.0)
        elem_matrix[:, 6] = range(start_index, start_index +(num_elem-1)*9, step=9.0)
    
        start_index = maximum(node_matrix[:,5]) + 2
        elem_matrix[:, 7] = range(start_index, start_index +(num_elem-1)*9, step=9.0)
    
        elem_matrix[:, 9] = ones(Float64, num_elem)
    
        start_index = maximum(node_matrix[:,5]) + 3
        elem_matrix[:, 11] = range(start_index, start_index +(num_elem-1)*9, step=9.0)
    
        start_index = maximum(node_matrix[:,5]) + 6
        elem_matrix[:, 14] = range(start_index, start_index +(num_elem-1)*9, step=9.0)
    
        start_index = maximum(node_matrix[:,5]) + 8
        elem_matrix[:, 15] = range(start_index, start_index +(num_elem-1)*9, step=9.0)
    
        start_index = maximum(node_matrix[:,5]) + 7
        elem_matrix[:, 18] = range(start_index, start_index +(num_elem-1)*9, step=9.0)
    
        start_index = maximum(node_matrix[:,5]) + 9
        elem_matrix[:, 19] = range(start_index, start_index +(num_elem-1)*9, step=9.0)
    
        elem_matrix[:, 20:27] .= Inf
    
        start_index = maximum(node_matrix[:,5]) + 4
        elem_matrix[:, 28] = range(start_index, start_index +(num_elem-1)*9, step=9.0)
    
        start_index = maximum(node_matrix[:,5]) + 5
        elem_matrix[:, 29] = range(start_index, start_index +(num_elem-1)*9, step=9.0)
    
        mastan_model["elem_info"] = elem_matrix
    
        uniload_matrix = zeros(Float64, num_elem, 16)
        uniload_matrix[:, 7:12] .= NaN
        uniload_matrix[:, 13] .= 6.50000000000000e-06
    
        mastan_model["uniload_info"] = uniload_matrix
    
        matwrite(joinpath(filepath, filename), mastan_model)
    
    end



end #module