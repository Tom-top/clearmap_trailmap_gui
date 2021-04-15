import os

import shutil
import numpy as np

import Parameters
import ClearMap.Alignment.Resampling as rs
import ClearMap.Alignment.Elastix as elx
import ClearMap.IO as io
from ClearMap.Alignment.Resampling import resamplePoints
from ClearMap.Alignment.Elastix import transformPoints
from ClearMap.Analysis.Voxelization import voxelize

from TrailMap.inference import segment_brain
from TrailMap.models import model

########################################################################################################################
# Setting generic file paths
########################################################################################################################

working_directory = "/raid/thomas.topilko/Grace_Trailmap_Analysis/trailmap_test_2"
auto_data = os.path.join(working_directory, "auto/autofluo.tif")
projections_data = os.path.join(working_directory, "projections/\d{4}.tif")

clearmap_gui_path = "/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek"
clearmap_path = os.path.join(clearmap_gui_path, "ClearMap")
clearmap_ressources_path = os.path.join(clearmap_gui_path, "ClearMap_Ressources")
template_file_path = os.path.join(clearmap_ressources_path, "25um_Autofluo_Reference/template_horizontal_1-289_25.tif")
trailmap_path = os.path.join(clearmap_gui_path, "TrailMap")
trailmap_weights_path = os.path.join(trailmap_path, "data/model-weights/trailmap_model.hdf5")
trailmap_result_directory = os.path.join(working_directory, "trailmap_results")
if not os.path.exists(trailmap_result_directory):
    os.mkdir(trailmap_result_directory)
trailmap_segmentation_result_directory = os.path.join(trailmap_result_directory, "segmentation")
if os.path.exists(trailmap_segmentation_result_directory):
    print(trailmap_segmentation_result_directory + " already exists. Will be overwritten")
    shutil.rmtree(trailmap_segmentation_result_directory)

########################################################################################################################
# Resample Auto and Projections data
########################################################################################################################

resampling_parameter_auto = Parameters.ResamplingParameterAutoFluo.copy()
resampling_parameter_auto["source"] = auto_data # The raw projections data
resampling_parameter_auto["sink"] = os.path.join(working_directory, "autofluo_resampled.tif")
resampling_parameter_auto["resolutionSource"] = (5, 5, 5)

rs.resampleData(**resampling_parameter_auto)

resampling_parameter_projections = Parameters.ResamplingParametercFos.copy()
resampling_parameter_projections["source"] = projections_data # The raw auto data
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
template_to_auto_directory = elx.alignData(**template_alignment_parameter)

########################################################################################################################
# TrailMap segmentation
########################################################################################################################

trailmap_model = model.get_net()
trailmap_model.load_weights(trailmap_weights_path)
segment_brain(os.path.basedir(projections_data), trailmap_segmentation_result_directory, trailmap_model)

########################################################################################################################
# Thresholding TrailMap segmentation result
########################################################################################################################

raw_trailmap_segmentation = io.readData(os.path.join(trailmap_segmentation_result_directory, "\d{4}.tif"))

# Thresholding
threshold = 0.7 # Threshold parameter below which all values in the TrailMap segmentation will be discarded
projections_mask = (raw_trailmap_segmentation < threshold)
projections_masked = raw_trailmap_segmentation.copy()
projections_masked[projections_mask] = 0
io.writeData(os.path.join(trailmap_result_directory, "projections_masked.tif"), projections_masked)

# Binarizing
projections_bin = projections_masked.copy()
projections_bin[~projections_mask] = 1
io.writeData(os.path.join(trailmap_result_directory, "projections_bin.tif"), projections_bin)

########################################################################################################################
# Calculating and saving coordinates of all positive points in the binarized TrailMap result
########################################################################################################################

point_coordinates = np.argwhere(projections_bin == 1)
points_coordinates_file = os.path.join(trailmap_result_directory, "point_coordinates.npy")
np.save(points_coordinates_file, point_coordinates)

params_resample_points = Parameters.ResamplePointsParameters.copy()
params_resample_points['pointSource'] = points_coordinates_file
params_resample_points['dataSizeSource'] = projections_file
points = resamplePoints(**params_resample_points)
points = transformPoints(points, transformDirectory=projections_to_auto_directory,
                         indices=False, resultDirectory=None, binary=False)
points = transformPoints(points, transformDirectory=template_to_auto_directory,
                         indices=False, resultDirectory=None, binary=False)
io.writePoints(os.path.join(trailmap_result_directory, "transformed_point_coordinates.npy"), points)

########################################################################################################################
# Creating the heatmap from point coordinates
########################################################################################################################

vox = voxelize(points, template_file_path, **Parameters.voxelizeParameter)
io.writeData(os.path.join(trailmap_result_directory, "heatmap_projections.tif"), vox.astype('int32'))