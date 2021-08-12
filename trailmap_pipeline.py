import os
import sys

import numpy as np
import pandas as pd
from shutil import copyfile, rmtree
import fileinput

import Parameters
import ClearMap.Alignment.Resampling as rs
import ClearMap.Alignment.Elastix as elx
import ClearMap.IO as io
from ClearMap.Alignment.Resampling import resamplePoints
from ClearMap.Alignment.Elastix import transformPoints
from ClearMap.Analysis.Voxelization import voxelize
import ClearMap.Tomek_Utilities as utils

from TrailMap.inference import segment_brain
from TrailMap.models import model

########################################################################################################################
# Setting generic file paths
########################################################################################################################

# samples_to_analyze = ["200117-340", "200117-344", "200117-373", "200117-399", "200117-404", "200117-407",
#                       "200117-412hydroceph", "200117-405", "200117-410", "200117-425", "200117-427", "200117-428"]
samples_to_analyze = ["200226-418", "200226-419", "200226-426", "200226-433"]

main_directory = "/raid/thomas.topilko/Grace_Trailmap_Analysis/grace_injection_data_barrel_to_barrel"
utils.check_all_paths_exist(main_directory, samples_to_analyze)

for s in samples_to_analyze:
    s = "200226-433"
    working_directory = os.path.join(main_directory, "{}".format(s))
    auto_data = os.path.join(working_directory, "auto/autofluo.tif")
    projections_data = os.path.join(working_directory, "projections/projections.tif")

    clearmap_gui_path = "/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek"
    clearmap_path = os.path.join(clearmap_gui_path, "ClearMap")
    clearmap_ressources_path = os.path.join(clearmap_gui_path, "ClearMap_Ressources")
    template_file_path = os.path.join(clearmap_ressources_path, "25um_Autofluo_Reference/template_horizontal_25.tif")
    trailmap_path = os.path.join(clearmap_gui_path, "TrailMap")
    trailmap_weights_path = os.path.join(trailmap_path, "data/model-weights/trailmap_model.hdf5")
    trailmap_result_directory = os.path.join(working_directory, "trailmap_results")
    if not os.path.exists(trailmap_result_directory):
        os.mkdir(trailmap_result_directory)
    trailmap_segmentation_result_directory = os.path.join(trailmap_result_directory, "segmentation")
    # if os.path.exists(trailmap_segmentation_result_directory):
    #     print(trailmap_segmentation_result_directory + " already exists. Will be overwritten")
    #     rmtree(trailmap_segmentation_result_directory)
    if not os.path.exists(trailmap_segmentation_result_directory):
        os.mkdir(trailmap_segmentation_result_directory)

    ########################################################################################################################
    # Resample Auto and Projections data
    ########################################################################################################################

    resampling_parameter_auto = Parameters.ResamplingParameterAutoFluo.copy()
    resampling_parameter_auto["source"] = auto_data # The raw auto data
    resampling_parameter_auto["sink"] = os.path.join(working_directory, "autofluo_resampled.tif")
    resampling_parameter_auto["resolutionSource"] = (5, 5, 5)

    rs.resampleData(**resampling_parameter_auto)

    resampling_parameter_projections = Parameters.ResamplingParametercFos.copy()
    resampling_parameter_projections["source"] = projections_data # The raw projections data
    resampling_parameter_projections["sink"] = os.path.join(working_directory, "projections_resampled.tif")
    resampling_parameter_projections["resolutionSource"] = (5, 5, 5)

    rs.resampleData(**resampling_parameter_projections)

    ########################################################################################################################
    # Align projections and auto
    ########################################################################################################################

    channels_alignment_parameter = Parameters.ChannelsAlignmentParameter.copy()
    channels_alignment_parameter["fixedImage"] = resampling_parameter_projections["sink"]
    channels_alignment_parameter["movingImage"] = resampling_parameter_auto["sink"]
    channels_alignment_parameter["affineParameterFile"] = os.path.join(clearmap_ressources_path, "Parameter_files/Par0000affine.txt")
    channels_alignment_parameter["resultDirectory"] = os.path.join(working_directory, "elastix_projections_to_auto")

    projections_to_auto_directory = os.path.join(working_directory, "elastix_projections_to_auto")
    projections_to_auto_directory = elx.alignData(**channels_alignment_parameter)

    ########################################################################################################################
    # Align template and auto
    ########################################################################################################################

    template_alignment_parameter = Parameters.TemplateAlignmentParameter.copy()
    template_alignment_parameter["fixedImage"] = resampling_parameter_auto["sink"]
    template_alignment_parameter["movingImage"] = template_file_path
    template_alignment_parameter["affineParameterFile"] = os.path.join(clearmap_ressources_path, "Parameter_files/Par0000affine_acquisition.txt")
    template_alignment_parameter["bSplineParameterFile"] = os.path.join(clearmap_ressources_path, "Parameter_files/Par0000bspline.txt")
    template_alignment_parameter["resultDirectory"] = os.path.join(working_directory, "elastix_template_to_auto")

    template_to_auto_directory = os.path.join(working_directory, "elastix_template_to_auto")
    template_to_auto_directory = elx.alignData(**template_alignment_parameter)

    ########################################################################################################################
    # Projection file into single files
    ########################################################################################################################

    single_file_projection_folder = os.path.join(working_directory, "projections/single_file_projections")
    if not os.path.exists(single_file_projection_folder):
        os.mkdir(single_file_projection_folder)
    single_file_projection_files = os.path.join(single_file_projection_folder, "\d{4}.tif")

    def stack_to_single_files(tiff_stack_path, sink_path):
        data = io.readData(tiff_stack_path)
        data = np.swapaxes(data, 0, -1)
        ndim = len(data.shape)
        n_planes = data.shape[1]
        for p in range(n_planes):
            if ndim == 3:
                io.writeData(os.path.join(sink_path, "{}.tif".format(str(p).zfill(4))),
                             data[:, p, :])
            elif ndim == 4:
                io.writeData(os.path.join(sink_path, "{}.tif".format(str(p).zfill(4))),
                             data[:, p, 2, :])

    stack_to_single_files(projections_data, single_file_projection_folder)
    ########################################################################################################################
    # TrailMap segmentation
    ########################################################################################################################

    trailmap_model = model.get_net()
    trailmap_model.load_weights(trailmap_weights_path)
    segment_brain.segment_brain(os.path.dirname(single_file_projection_files), trailmap_segmentation_result_directory, trailmap_model)

    ########################################################################################################################
    # Thresholding TrailMap segmentation result
    ########################################################################################################################

    raw_trailmap_segmentation = io.readData(os.path.join(trailmap_segmentation_result_directory, "seg-\d{4}.tif"))

    # Thresholding
    threshold = 0.7 # Threshold parameter below which all values in the TrailMap segmentation will be discarded
    projections_mask = (raw_trailmap_segmentation < threshold)
    projections_masked = raw_trailmap_segmentation.copy()
    projections_masked[projections_mask] = 0
    print("finished thresholding trailmap data for animal {}".format(s))
    # io.writeData(os.path.join(trailmap_result_directory, "projections_masked.tif"), projections_masked)

    # Binarizing
    projections_bin = projections_masked.copy()
    projections_bin[~projections_mask] = 1
    print("finished masking trailmap data for animal {}".format(s))
    io.writeData(os.path.join(trailmap_result_directory, "projections_bin.tif"), projections_bin)

    ########################################################################################################################
    # Calculating and saving coordinates of all positive points in the binarized TrailMap result
    ########################################################################################################################

    projections_bin = io.readData(os.path.join(trailmap_result_directory, "projections_bin.tif"))
    point_coordinates = np.argwhere(projections_bin == 1)
    points_coordinates_file = os.path.join(trailmap_result_directory, "point_coordinates.npy")
    print("finished computing coordinates for trailmap data for animal {}".format(s))
    np.save(points_coordinates_file, point_coordinates)

    import matplotlib.pyplot as plt
    import matplotlib.patches as patch
    size = io.dataSize(projections_bin)
    points_test = io.readPoints(points_coordinates_file)
    fig = plt.figure(dpi=300)
    ax0 = plt.subplot(111, aspect="equal")
    ax0.scatter(points_test[:,0][::100], points_test[:,1][::100], s=0.5, alpha=0.2)
    rect = patch.Rectangle((0, 0), size[0], size[1], linewidth=1, edgecolor='r', facecolor='none')
    ax0.add_patch(rect)
    plt.gca().invert_yaxis()
    plt.show()
    plt.gcf()

    # np.load(points_coordinates_file)
    params_resample_points = Parameters.ResamplePointsParameters.copy()
    params_resample_points['pointSource'] = points_coordinates_file
    params_resample_points['dataSizeSource'] = projections_data
    points = resamplePoints(**params_resample_points)

    # import matplotlib.pyplot as plt
    # size = io.dataSize(os.path.join(working_directory, "projections_resampled.tif"))
    # fig = plt.figure(dpi=300)
    # ax0 = plt.subplot(111, aspect="equal")
    # ax0.scatter(points[:, 0][::100], points[:, 1][::100], s=0.5, alpha=0.2)
    # rect = patch.Rectangle((0, 0), size[0], size[1], linewidth=1, edgecolor='r', facecolor='none')
    # ax0.add_patch(rect)
    # plt.gca().invert_yaxis()
    # plt.show()

    # transform_trailmap_to_auto = os.path.join(working_directory, "trailmap_results/transformix_trailmap_to_auto")
    points_trailmap_to_auto = transformPoints(points, transformDirectory=projections_to_auto_directory,
                                              indices=False, resultDirectory=None, binary=False)

    # import matplotlib.pyplot as plt
    # from matplotlib.widgets import MultiCursor
    # size = io.dataSize(os.path.join(working_directory, "autofluo_resampled.tif"))
    # raw_data = io.readData(os.path.join(working_directory, "autofluo_resampled.tif"))
    # orientation = [0, 2]
    # z_project = np.max(raw_data, axis=1)
    # # z_project = raw_data[:, 350, :]
    # fig = plt.figure(dpi=300)
    # ax0 = plt.subplot(1, 2, 1, aspect="equal")
    # ax0.scatter(points_trailmap_to_auto[:, orientation[0]][::100],
    #             points_trailmap_to_auto[:, orientation[1]][::100], s=0.5, alpha=0.2)
    # # rect = patch.Rectangle((0, 0), size[orientation[0]], size[orientation[1]], linewidth=1, edgecolor='r', facecolor='none')
    # # ax0.add_patch(rect)
    # ax0.set_xticks([])
    # ax0.set_yticks([])
    # plt.gca().invert_yaxis()
    # ax1 = plt.subplot(1, 2, 2, sharex=ax0, sharey=ax0)
    # ax1.imshow(np.swapaxes(z_project, 1, 0), vmin=0, vmax=10000)
    # ax1.set_xticks([])
    # ax1.set_yticks([])
    # multi = MultiCursor(fig.canvas, (ax0, ax1), horizOn=True, vertOn=True, color='r', lw=0.2)
    # plt.show()

    # sink = np.zeros(raw_data.shape)
    # for i in points_trailmap_to_auto:
    #     x, y, z = np.round(i)
    #     try:
    #         sink[int(x), int(y), int(z)] = 1
    #     except:
    #         print("{} out of bounds!".format([x, y, z]))
    # io.writeData(os.path.join(working_directory, "transformed_points_to_auto.tif"), sink)

    # transform_trailmap_to_template = os.path.join(working_directory, "trailmap_results/transformix_trailmap_to_template")
    points_trailmap_to_template = transformPoints(points_trailmap_to_auto, transformDirectory=template_to_auto_directory,
                                                  indices=False, resultDirectory=None, binary=False)

    # import matplotlib.pyplot as plt
    # size = io.dataSize(template_file_path)
    # fig = plt.figure(dpi=300)
    # ax0 = plt.subplot(111, aspect="equal")
    # orientation = [0, 1]
    # ax0.scatter(points_trailmap_to_template[:, orientation[0]][::100],
    #             points_trailmap_to_template[:, orientation[1]][::100], s=0.5, alpha=0.2)
    # rect = patch.Rectangle((0, 0), size[orientation[0]], size[orientation[1]], linewidth=1, edgecolor='r', facecolor='none')
    # ax0.add_patch(rect)
    # plt.gca().invert_yaxis()
    # plt.show()

    print("finished transforming points for trailmap data for animal {}".format(s))
    io.writePoints(os.path.join(trailmap_result_directory, "transformed_point_coordinates.npy"), points_trailmap_to_template)

    ########################################################################################################################
    # Creating the heatmap from point coordinates
    ########################################################################################################################
# for s in samples_to_analyze:
#     working_directory = os.path.join(main_directory, "{}".format(s))
#     trailmap_result_directory = os.path.join(working_directory, "trailmap_results")
#     points = io.readPoints(os.path.join(trailmap_result_directory, "transformed_point_coordinates.npy"))
#     clearmap_gui_path = "/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek"
#     clearmap_ressources_path = os.path.join(clearmap_gui_path, "ClearMap_Ressources")
#     template_file_path = os.path.join(clearmap_ressources_path, "25um_Autofluo_Reference/template_horizontal_1-289_25.tif")
    points_trailmap_to_template = io.readPoints(os.path.join(trailmap_result_directory, "transformed_point_coordinates.npy"))
    voxelization_parameters = Parameters.voxelizeParameter.copy()
    # voxelization_parameters["size"] = (20, 20, 20)
    voxelization_parameters["size"] = (15, 15, 15)
    vox = voxelize(points_trailmap_to_template, template_file_path, **voxelization_parameters)

    import matplotlib.pyplot as plt
    fig = plt.figure(dpi=300)
    plt.scatter(vox[:, 0][::100], vox[:, 2][::100], s=0.5, alpha=0.2)
    plt.gca().invert_yaxis()
    plt.show()

    print("finished voxelizing trailmap data for animal {}".format(s))
    io.writeData(os.path.join(trailmap_result_directory, "heatmap_projections.tif"), vox.astype('int32'))


########################################################################################################################
########################################################################################################################
# Coloring TrailMap result
########################################################################################################################
########################################################################################################################

working_directory = '/raid/thomas.topilko/Grace_Trailmap_Analysis/grace_injection_data/200117-336'
annotation_file_path = os.path.join(clearmap_ressources_path, "Regions_annotations/ABA_25um_annotation_horizontal.tif")
# annotation_file_path = "/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek/ClearMap_Ressources/25um_Autofluo_Reference/template_horizontal_25.tif"
trailmap_coloring_directory = os.path.join(working_directory, "trailmap_coloring")
if not os.path.exists(trailmap_coloring_directory):
    os.mkdir(trailmap_coloring_directory)

# channels_alignment_parameter_reverse = Parameters.ChannelsAlignmentParameter.copy()
# channels_alignment_parameter_reverse["fixedImage"] = os.path.join(working_directory, "autofluo_resampled.tif")
# channels_alignment_parameter_reverse["movingImage"] = os.path.join(working_directory, "projections_resampled.tif")
# channels_alignment_parameter_reverse["affineParameterFile"] = os.path.join(clearmap_ressources_path, "Parameter_files/Par0000affine.txt")
# channels_alignment_parameter_reverse["resultDirectory"] = os.path.join(trailmap_coloring_directory, "elastix_auto_to_projections")
#
# projections_to_auto_directory = elx.alignData(**channels_alignment_parameter_reverse)
#
# template_alignment_parameter_reverse = Parameters.TemplateAlignmentParameter.copy()
# template_alignment_parameter_reverse["fixedImage"] = template_file_path
# template_alignment_parameter_reverse["movingImage"] = os.path.join(working_directory, "autofluo_resampled.tif")
# template_alignment_parameter_reverse["affineParameterFile"] = os.path.join(clearmap_ressources_path, "Parameter_files/Par0000affine.txt")
# template_alignment_parameter_reverse["bSplineParameterFile"] = os.path.join(clearmap_ressources_path, "Parameter_files/Par0000bspline.txt")
# template_alignment_parameter_reverse["resultDirectory"] = os.path.join(trailmap_coloring_directory, "elastix_auto_to_template")
#
# projections_to_auto_directory = elx.alignData(**template_alignment_parameter_reverse)

copyfile(os.path.join(template_alignment_parameter["resultDirectory"],'TransformParameters.1.txt'),
         os.path.join(template_alignment_parameter["resultDirectory"],'TransformParameters_no_interpolation.1.txt'))
transform_file = os.path.join(template_alignment_parameter["resultDirectory"], 'TransformParameters_no_interpolation.1.txt')
for line in fileinput.input([transform_file], inplace=True):
   if 'InitialTransformParametersFileName' in line:
       line = line.replace(line, '(InitialTransformParametersFileName "{}")'.
                           format(os.path.join(template_alignment_parameter["resultDirectory"],
                                               'TransformParameters.0.txt')))
   line = line.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)')
   sys.stdout.write(line)

transform_annotation_to_auto_dir = os.path.join(trailmap_coloring_directory, "transformix_auto_to_annotation")
transform_annotation_to_auto = os.path.join(transform_annotation_to_auto_dir, "transformed_auto_to_annotation.tif")
elx.transformData(annotation_file_path, sink=transform_annotation_to_auto,
                  transformParameterFile=transform_file, resultDirectory=transform_annotation_to_auto_dir)

copyfile(os.path.join(channels_alignment_parameter["resultDirectory"],'TransformParameters.0.txt'),
         os.path.join(channels_alignment_parameter["resultDirectory"],'TransformParameters_no_interpolation.0.txt'))
transform_file = os.path.join(channels_alignment_parameter["resultDirectory"], 'TransformParameters_no_interpolation.0.txt')
for line in fileinput.input([transform_file], inplace=True):
   line = line.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)')
   sys.stdout.write(line)

transform_annotation_to_projections_dir = os.path.join(trailmap_coloring_directory, "transformix_projections_to_annotation")
transform_annotation_to_projections = os.path.join(transform_annotation_to_projections_dir, "transformed_projections_to_annotation.tif")
elx.transformData(transform_annotation_to_auto, sink=transform_annotation_to_projections,
                  transformParameterFile=transform_file, resultDirectory=transform_annotation_to_projections_dir)

resampling_parameter_projections_10m = Parameters.ResamplingParametercFos.copy()
resampling_parameter_projections_10m["source"] = os.path.join(trailmap_result_directory, "projections_bin.tif")
resampling_parameter_projections_10m["sink"] = os.path.join(trailmap_coloring_directory, "projections_resampled_10m.tif")
resampling_parameter_projections_10m["resolutionSource"] = (5, 5, 5)
resampling_parameter_projections_10m["resolutionSink"] = (10, 10, 10)

rs.resampleData(**resampling_parameter_projections_10m)

resampling_parameter_transformed_annotation_10m = Parameters.ResamplingParametercFos.copy()
resampling_parameter_transformed_annotation_10m["source"] = transform_annotation_to_projections
resampling_parameter_transformed_annotation_10m["sink"] = os.path.join(trailmap_coloring_directory, "transformed_template_resampled_10m.tif")
resampling_parameter_transformed_annotation_10m["resolutionSource"] = (25, 25, 25)
resampling_parameter_transformed_annotation_10m["resolutionSink"] = (10, 10, 10)
resampling_parameter_transformed_annotation_10m["interpolation"] = None

rs.resampleData(**resampling_parameter_transformed_annotation_10m)

samples_to_analyze = ["200226-433"]
for s in samples_to_analyze:
    working_directory = os.path.join(main_directory, "{}".format(s))
    auto_data = os.path.join(working_directory, "auto/autofluo.tif")
    projections_data = os.path.join(working_directory, "projections/projections.tif")

    clearmap_gui_path = "/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek"
    clearmap_path = os.path.join(clearmap_gui_path, "ClearMap")
    clearmap_ressources_path = os.path.join(clearmap_gui_path, "ClearMap_Ressources")
    template_file_path = os.path.join(clearmap_ressources_path, "25um_Autofluo_Reference/template_horizontal_25.tif")
    trailmap_path = os.path.join(clearmap_gui_path, "TrailMap")
    trailmap_weights_path = os.path.join(trailmap_path, "data/model-weights/trailmap_model.hdf5")
    trailmap_result_directory = os.path.join(working_directory, "trailmap_results")
    if not os.path.exists(trailmap_result_directory):
        os.mkdir(trailmap_result_directory)
    trailmap_segmentation_result_directory = os.path.join(trailmap_result_directory, "segmentation")
    # if os.path.exists(trailmap_segmentation_result_directory):
    #     print(trailmap_segmentation_result_directory + " already exists. Will be overwritten")
    #     rmtree(trailmap_segmentation_result_directory)
    if not os.path.exists(trailmap_segmentation_result_directory):
        os.mkdir(trailmap_segmentation_result_directory)

    resampling_parameter_transformed_annotation_10m = Parameters.ResamplingParametercFos.copy()
    resampling_parameter_transformed_annotation_10m["source"] = os.path.join(trailmap_result_directory, "projections_bin.tif")
    resampling_parameter_transformed_annotation_10m["sink"] = os.path.join(trailmap_result_directory, "projections_bin_resampled_10m.tif")
    resampling_parameter_transformed_annotation_10m["resolutionSource"] = (5, 5, 5)
    resampling_parameter_transformed_annotation_10m["resolutionSink"] = (25, 25, 25)
    resampling_parameter_transformed_annotation_10m["interpolation"] = None

    rs.resampleData(**resampling_parameter_transformed_annotation_10m)



regions_ids_path = os.path.join(clearmap_ressources_path, "Regions_annotations/regions_IDs.csv")
regions_ids = pd.read_csv(regions_ids_path)
projections_10m = io.readData(resampling_parameter_projections_10m["sink"])
projections_10m_rgb = np.stack((projections_10m,)*3, axis=-1)
mask_projections_10m = projections_10m == 1
annotation_10m = io.readData(resampling_parameter_transformed_annotation_10m["sink"])
annotation_10m = annotation_10m[:, :, :-1]
gray_values = np.unique(annotation_10m)

for gray_value in gray_values:
    if gray_value > 0:
        print("{}/{}".format(gray_value, gray_values[-1]))
        # positions_in_annotation = np.argwhere(annotation_10m == gray_value)
        positions_in_annotation = annotation_10m == gray_value
        pos_in_regions_ids = regions_ids.index[regions_ids['id'] == gray_value].tolist()
        if not len(pos_in_regions_ids) == 1:
            print("id: {} has multiple entries".format(gray_value))
        pos_in_regions_ids = pos_in_regions_ids[0]
        red_color_value = regions_ids['red'][pos_in_regions_ids]
        green_color_value = regions_ids['green'][pos_in_regions_ids]
        blue_color_value = regions_ids['blue'][pos_in_regions_ids]
        color = (red_color_value, green_color_value, blue_color_value)
        mask = np.logical_and(positions_in_annotation, mask_projections_10m)
        projections_10m_rgb[mask] = color

io.writeData(os.path.join(trailmap_coloring_directory, "colored_projections_10m.tif"), projections_10m_rgb)

########################################################################################################################
########################################################################################################################
# Generate cortex atlas
########################################################################################################################
########################################################################################################################

full_annotation_file_path = "/raid/thomas.topilko/Grace_Trailmap_Analysis/grace_injection_data/200117-340/trailmap_coloring/transformed_annotation_resampled_10m.tif"
regions_ids_cortex_file_path = "/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek/ClearMap_Ressources/Regions_annotations/regions_IDs_cortex.xlsx"

annotation_data = io.readData(full_annotation_file_path)
regions_ids_data = pd.read_excel(regions_ids_cortex_file_path)
ids_to_keep = regions_ids_data["id"]

mask_cortex = np.isin(annotation_data, ids_to_keep)
cortex_annotation_data = np.where(mask_cortex, annotation_data, 0)
io.writeData(os.path.join("/raid/thomas.topilko/Grace_Trailmap_Analysis/grace_injection_data/200117-340/trailmap_coloring", "transformed_annotation_resampled_10m_cortex.tif"),
             cortex_annotation_data)






















