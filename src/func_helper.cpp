#include "func_helper.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/bounding_box.h>
#include <igl/signed_distance.h>
using namespace Eigen;
AlignedBox3d generate_bounding_box(const MatrixXd& vertices_matrix);
MatrixXd scale_and_transform_vertices(const MatrixXd& vertices_matrix, const AlignedBox3d& desired_bounding_box, bool is_output_same_size_as_original);
MatrixXd mesh_axises_align(const MatrixXd& vertices_matrix);
void grid_from_mesh(const MatrixXd& vertices_matrix, MatrixXd& voxels_centers, Vector3i& grid_dimensions_in_voxel, Vector3d& voxel_size, MatrixXd& center_matrix);
AlignedBox3d generate_alinged_box(const MatrixXd& vertices_matrix);


Vector3d compute_diagonal_of_box(Vector3d first_corner, Vector3d second_corner)
{
	return first_corner - second_corner;
}

MatrixXd scale_and_transform_vertices(const MatrixXd& vertices_matrix, const AlignedBox3d& desired_bounding_box, bool is_output_same_size_as_original)
{
	MatrixXd result;
	int resize_factor = vertices_matrix.rows();
	result.resize(resize_factor, 3);
	//calculate the diagonals (bottom left floor - top right ceil) of the input and the output grids
	Vector3d bottom_left_floor_corner_target_box;
	bottom_left_floor_corner_target_box = desired_bounding_box.corner(desired_bounding_box.BottomLeftFloor);
	Vector3d top_right_ceil_corner_target_box;
	top_right_ceil_corner_target_box = desired_bounding_box.corner(desired_bounding_box.TopRightCeil);
	AlignedBox<double, 3> input_bounding_box;
	input_bounding_box = generate_alinged_box(vertices_matrix);
	Vector3d bottom_left_floor_corner_source_box;
	bottom_left_floor_corner_source_box = input_bounding_box.corner(input_bounding_box.BottomLeftFloor);
	Vector3d top_right_ceil_corner_source_box;
	top_right_ceil_corner_source_box = input_bounding_box.corner(input_bounding_box.TopRightCeil);
	Vector3d diagonal_of_input_box;
	diagonal_of_input_box = compute_diagonal_of_box(bottom_left_floor_corner_source_box, top_right_ceil_corner_source_box);
	Vector3d diagonal_of_output_box;
	diagonal_of_output_box = compute_diagonal_of_box(bottom_left_floor_corner_target_box, top_right_ceil_corner_target_box);
	// resize vertivces
	Matrix3d scaling_factor;
	scaling_factor = Matrix3d::Identity();
	if (false == is_output_same_size_as_original) scaling_factor.diagonal() << diagonal_of_output_box.cwiseQuotient(diagonal_of_input_box).array();
	result = vertices_matrix* scaling_factor;
	diagonal_of_input_box = scaling_factor * bottom_left_floor_corner_source_box - scaling_factor * top_right_ceil_corner_source_box;//perform rotation
	Matrix3d rotation_matrix;
	rotation_matrix = igl::rotation_matrix_from_directions(diagonal_of_input_box.normalized(), diagonal_of_output_box.normalized());
	result = result * rotation_matrix.transpose();
	Vector3d translation_vector;
	translation_vector = bottom_left_floor_corner_target_box - rotation_matrix * scaling_factor * bottom_left_floor_corner_source_box; 
	for (int i = 0; i < result.rows(); i++) 
		result.row(i) += translation_vector.transpose();
	return result;
}
/// <summary>
/// generating the bounding box, given the vertices
/// </summary>
/// <param name="vertices_matrix"></param>
/// <returns></returns>
AlignedBox3d generate_bounding_box(const MatrixXd& vertices_matrix)
{
	MatrixXd vertices_helper;
	MatrixXd faces_helper;
	igl::bounding_box(vertices_matrix, vertices_helper, faces_helper);
	AlignedBox3d res;
	for (int i = 0; i < vertices_matrix.rows(); i++)
		res.extend(vertices_matrix.row(i).transpose());
	return res;
}
MatrixXd mesh_axises_align(const MatrixXd& vertices_matrix)
{
	MatrixXd res, aligned_grid;
	AlignedBox3d original_box, aligned_box;
	Vector3d dimension_vector;
	original_box = generate_alinged_box(vertices_matrix);
	dimension_vector = Vector3d::Constant(10.0); //sets all coefficients to 10
	igl::grid(dimension_vector, aligned_grid); ;//creating the grid
	aligned_box = generate_alinged_box(aligned_grid);
	res = scale_and_transform_vertices(vertices_matrix, aligned_box, true); 
	return res;
}
void grid_from_mesh(const MatrixXd& vertices_matrix, MatrixXd& voxels_centers, Vector3i& grid_dimensions_in_voxel, Vector3d& voxel_size, MatrixXd& center_matrix)
{
	Vector3d scaling_factor;
	Vector3i ones_vector(1, 1, 1);
	AlignedBox3d bounding_box;
	bounding_box = generate_bounding_box(vertices_matrix);
	RowVector3i grid_dimension_vector;
	igl::voxel_grid(bounding_box, 55, 0, center_matrix, grid_dimension_vector); //Creating the voxel grid using igl. I bound 55 voxels at each side of the grid.
	grid_dimensions_in_voxel = grid_dimension_vector.transpose();
	AlignedBox3d box_of_centers;
	box_of_centers = generate_alinged_box(center_matrix);
	voxel_size = box_of_centers.sizes().cwiseQuotient((grid_dimensions_in_voxel - ones_vector).cast<double>());
	MatrixXd grid;
	Vector3d corners_dimensions = (grid_dimensions_in_voxel + ones_vector).cast<double>();
	igl::grid(corners_dimensions, grid); // Construct vertices of a regular grid
	Vector3d bottom_left_floor;
	scaling_factor = voxel_size / 2.0;
	bottom_left_floor = box_of_centers.corner(box_of_centers.BottomLeftFloor) - (scaling_factor); //changing the original grid into the grid of voxels
	Vector3d top_right_ceil;
	top_right_ceil = box_of_centers.corner(box_of_centers.TopRightCeil) + (scaling_factor);
	AlignedBox3d final_box(bottom_left_floor, top_right_ceil);
	voxels_centers = scale_and_transform_vertices(grid, final_box, false);
}

AlignedBox3d generate_alinged_box(const MatrixXd& vertices_matrix)
{
	AlignedBox3d res_box;
	for (int i = 0; i < vertices_matrix.rows(); i++) {
		res_box.extend(vertices_matrix.row(i).transpose());
	}
	return res_box;
}