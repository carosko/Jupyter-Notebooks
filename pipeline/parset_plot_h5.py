! working_directory              = /data009/scratch/carosko/work
! job_directory                  = /data009/scratch/carosko/job
! res_directory                  = /data009/scratch/carosko/res
#################################################
! plot_solH5_script              = /home/carosko/lbbw/comp_cha/script/jupy_l.py
#################################################
! target_input_path              = /data020/scratch/lb_bw/test_losoto_delay_calibration/noCS_circular/
#/data020/scratch/lb_bw/test_losoto_delay_calibration/16sec_1SB_2comp_uvcutCS/
! target_input_pattern           = *.h5

#! chess_switch                   = True
#! chess_parameter_a              = 1600.
#! chess_parameter_b              = 600.


pipeline.steps                  = [mapfiles,  combine_data_tar_map, plot_solH5]

# mapfile
mapfiles.control.kind            = plugin
mapfiles.control.type            = createMapfile
mapfiles.control.method          = mapfile_from_folder
mapfiles.control.mapfile_dir     = input.output.mapfile_dir
mapfiles.control.filename        = mapfiles.mapfile
mapfiles.control.folder          = {{ target_input_path }}
mapfiles.control.pattern         = {{ target_input_pattern }}

# combine all entries into one mapfile 
combine_data_tar_map.control.kind            =   plugin
combine_data_tar_map.control.type            =   createMapfile
combine_data_tar_map.control.method          =   mapfile_all_to_one
combine_data_tar_map.control.mapfile_dir     =   input.output.mapfile_dir
combine_data_tar_map.control.filename        =   combine_data_tar_map.mapfile
combine_data_tar_map.control.mapfile_in      =   mapfiles.output.mapfile




# plot_solH5
plot_solH5.control.type             = pythonplugin
plot_solH5.control.executable       = {{ plot_solH5_script }}
plot_solH5.argument.h5_list         = combine_data_tar_map.output.mapfile
plot_solH5.argument.output_dir      = {{ res_directory }}
#plot_solH5.argument.output_clipping = {{ res_directory }}/clipping.txt
#chessboard.argument.parameter_a = {{ chess_parameter_a }}
#chessboard.argument.parameter_b = {{ chess_parameter_b }}
#chessboard.argument.output_name = {{ working_directory }}/chessboard.txt
#chessboard.argument.switch      = {{ chess_switch }}
# chessboard.argument.cen       = hard-coded in the script
# chessboard.argument.throw     = hard-coded in the script






