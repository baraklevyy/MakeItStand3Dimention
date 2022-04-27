#include "carving.h"
#include <igl/voxel_grid.h>
#include "func_helper.h"
#include <limits>
#include <iostream>
#include <igl/centroid.h>
#include <igl/per_face_normals.h>
#include <algorithm>
#include <stdio.h>
#include <igl/bounding_box.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>
using namespace Eigen;
typedef struct axes_max_min
{
	double max_axis_x;
    double min_axis_x;
	double max_axis_y;
	double min_axis_y;
	double max_axis_z;
	double min_axis_z;
}axes_max_min;

typedef struct row_vecotr_container
{
	RowVector3d A;
	RowVector3d B;
	RowVector3i triangle;
	RowVector3d N;
	RowVector3d X2;
	RowVector3d X1;
	RowVector3d X0;
}row_vec_container;

typedef struct num_container
{
	double a;
	double b;
	double c;
	double s;
	double area;
}num_container;

const int ERROR_INTERSECT = 4;
std::vector<int> get_index_of_adjacent(double infinitesimal, double border_of_current, double max, double min, double voxel_size);
bool is_border_voxel(int x_index, int y_index, int z_index, const Vector3i vector);
bool in_range(std::vector<int>& neighabors_of_x, std::vector<int>& neighabors_of_y, std::vector<int>& neighabors_of_z, int outer_itration, int middle_iteration, int inner_iteration);
void corner_iteration(Vector3d& voxel_size, std::vector<int>& neighabors_of_x, std::vector<int>& neighabors_of_y, std::vector<int>& neighabors_of_z,
	List3d<Voxel>& interior_voxel_mesh, const double infinitesimal, axes_max_min* max_min, const MatrixXd voxel_corners);
bool necessary_to_append_face(const Voxel& first_voxel, const Voxel& second_voxel);
MatrixXi triangulate_matrix(std::vector<int> intersaction_vector);
void triangulte_common_side_voxels(MatrixXi& pushed_matrix, const Voxel& first_voxel, const Voxel& second_voxel, const MatrixXd& vertices_matrix);
void com_step(Vector3d& com_g, const Voxel& added_voxel, const Voxel& neig_voxel, const MatrixXd& vertices_matrix, double mass);
void center_of_mass(const MatrixXd& vertices_matrix, const MatrixXi& faces_matrix, Vector3d& centers_vec, double mass);


const int NUMBER_OF_NEIGHABORS = 2;


bool is_border_voxel(int x_index, int y_index, int z_index, const Vector3i vector)
{
	if ((x_index > 0 && x_index < vector(0) - 1) && (y_index > 0 && y_index < vector(1) - 1) && (z_index > 0 && z_index < vector(2) - 1)) return true;
	return false;
}
bool in_range(std::vector<int>& neighabors_of_x, std::vector<int>& neighabors_of_y, std::vector<int>& neighabors_of_z, int outer_itration, int middle_iteration, int inner_iteration)
{
	if (neighabors_of_x[outer_itration] != -1 && neighabors_of_y[middle_iteration] != -1 && neighabors_of_z[inner_iteration] != -1) return true;
	return false;
}
void corner_iteration(Vector3d& voxel_size, std::vector<int>& neighabors_of_x, std::vector<int>& neighabors_of_y, std::vector<int>& neighabors_of_z,
	List3d<Voxel>& interior_voxel_mesh, const double infinitesimal, axes_max_min* max_min, const MatrixXd voxel_corners)
{
	RowVector3d current_iteration_corner;
	for (int rowIdx = 0; rowIdx < voxel_corners.rows(); rowIdx++)
	{
		current_iteration_corner.setZero();
		current_iteration_corner = voxel_corners.row(rowIdx);
		neighabors_of_x = get_index_of_adjacent(infinitesimal, current_iteration_corner(0), max_min->max_axis_x, max_min->min_axis_x, voxel_size(0));
		neighabors_of_y = get_index_of_adjacent(infinitesimal, current_iteration_corner(1), max_min->max_axis_y, max_min->min_axis_y, voxel_size(1));
		neighabors_of_z = get_index_of_adjacent(infinitesimal, current_iteration_corner(2), max_min->max_axis_z, max_min->min_axis_z, voxel_size(2));

		for (int i = 0; i < NUMBER_OF_NEIGHABORS; i++) {
			for (int j = 0; j < NUMBER_OF_NEIGHABORS; j++) {
				for (int k = 0; k < NUMBER_OF_NEIGHABORS; k++)
					if (in_range(neighabors_of_x, neighabors_of_y, neighabors_of_z, i, j, k)) interior_voxel_mesh(neighabors_of_x[i], neighabors_of_y[j], neighabors_of_z[k]).voxel_edge_vertices_index.insert(rowIdx);
			}
		}
	}
}
std::vector<int> get_index_of_adjacent(double infinitesimal, double border_of_current, double max, double min, double voxel_size)
{
	std::vector<int> res(2);
	assert(border_of_current >= min - infinitesimal && border_of_current <= max + infinitesimal);
	int limit_index = std::round((border_of_current - min) / voxel_size);
	assert(std::abs(border_of_current - (min + (limit_index) * voxel_size)) < infinitesimal);
	int maximum_index = std::round((max - min) / voxel_size);
	if (limit_index == 0) res.assign({ 0, -1 });
	else if (limit_index == maximum_index) res.assign({ limit_index - 1, -1 });
	else res.assign({ limit_index - 1, limit_index});
	return res;
}

/// <summary>
/// constructor and initiator
/// </summary>
/// <param name="gravity"></param>
/// <param name="equilibrium_point"></param>
/// <param name="vertices_matrix"></param>
/// <param name="faces_matrix"></param>
/// <param name="plane_of_vertices"></param>
/// <param name="plane_of_faces"></param>
Carving_State::Carving_State(const Vector3d& gravity, const Vector3d& equilibrium_point, const MatrixXd & vertices_matrix, const MatrixXi& faces_matrix,
	const MatrixXd &plane_of_vertices, const MatrixXi &plane_of_faces)
{
	this->gravity_direction = gravity;
	this->equilibrium_point = equilibrium_point;
	this->vertices_plane = plane_of_vertices;
	this->faces_plane = plane_of_faces;
	this->optimal_com = equilibrium_point.transpose();
	this->interior_voxel_mesh = List3d<Voxel>();
	this->number_of_voxels_to_display = 0;
	this->number_of_iter = 0;
	double com_mass;
	igl::centroid(vertices_matrix, faces_matrix, this->current_com, com_mass);
	this->object_mass = (com_mass);
	MatrixXd matrix_of_centers;
	grid_from_mesh(vertices_matrix, this->voxel_corners, this->voxel_dimensions, this->voxel_size, matrix_of_centers);
	AlignedBox3d border_of_grid;
	border_of_grid = generate_alinged_box(voxel_corners);
	axes_max_min max_min;
	max_min.min_axis_x = voxel_corners.col(0).minCoeff();
	max_min.max_axis_x = voxel_corners.col(0).maxCoeff();
	max_min.min_axis_y = voxel_corners.col(1).minCoeff();
	max_min.max_axis_y = voxel_corners.col(1).maxCoeff();
	max_min.min_axis_z = voxel_corners.col(2).minCoeff();
	max_min.max_axis_z = voxel_corners.col(2).maxCoeff();
	VectorXd dist, min_dist_faces;
	MatrixXd min_dist_point, min_dist_norm;
	constexpr double minimum_limit = std::numeric_limits<double>::min();
	constexpr double maximum_limit = std::numeric_limits<double>::max();
	igl::signed_distance(matrix_of_centers, vertices_matrix, faces_matrix, igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER, minimum_limit, maximum_limit, dist, min_dist_faces, min_dist_point, min_dist_norm); //calculated the distance ti mesh in voxels
	interior_voxel_mesh.resize(this->voxel_dimensions(0), this->voxel_dimensions(1), this->voxel_dimensions(2)); //set the interior mesh dimensions
	RowVector3d center_of_curr_iteration;
	int index_x_axis, index_y_axis, index_z_axis;
	int num_of_rows = matrix_of_centers.rows();
	for (int i = 0; i < num_of_rows; i++)
	{
		center_of_curr_iteration = matrix_of_centers.row(i);
		assert(center_of_curr_iteration(0) <= max_min.max_axis_x && center_of_curr_iteration(0) >= max_min.min_axis_x);
		index_x_axis = std::floor((center_of_curr_iteration(0) - max_min.min_axis_x) / this->voxel_size(0));
		assert(center_of_curr_iteration(1) <= max_min.max_axis_y && center_of_curr_iteration(1) >= max_min.min_axis_y);
		index_y_axis = std::floor((center_of_curr_iteration(1) - max_min.min_axis_y) / this->voxel_size(1));
		assert(center_of_curr_iteration(2) <= max_min.max_axis_z && center_of_curr_iteration(2) >= max_min.min_axis_z);
		index_z_axis = std::floor((center_of_curr_iteration(2) - max_min.min_axis_z) / this->voxel_size(2));
		interior_voxel_mesh(index_x_axis, index_y_axis, index_z_axis).center = center_of_curr_iteration;
		interior_voxel_mesh(index_x_axis, index_y_axis, index_z_axis).index_in_grid_s.axis_x_index = index_x_axis;
		interior_voxel_mesh(index_x_axis, index_y_axis, index_z_axis).index_in_grid_s.axis_y_index = index_y_axis;
		interior_voxel_mesh(index_x_axis, index_y_axis, index_z_axis).index_in_grid_s.axis_z_index = index_z_axis;
		interior_voxel_mesh(index_x_axis, index_y_axis, index_z_axis).closet_dist_to_mesh = dist(i);
		if (dist(i) < 0.0 && std::abs(dist(i)) >= voxel_size(0))
		{
			bool border_voxel = is_border_voxel(index_x_axis, index_y_axis, index_z_axis, voxel_dimensions);
			if(border_voxel) interior_voxel_mesh(index_x_axis, index_y_axis, index_z_axis).is_interior_voxel = true; //mark as border voxel. I don't calculate them within the group pf voxels that can be carved
		}
	}
	std::vector<int> neighabors_of_x(2);
	std::vector<int> neighabors_of_y(2);
	std::vector<int> neighabors_of_z(2);
	corner_iteration(voxel_size, neighabors_of_x, neighabors_of_y, neighabors_of_z, this->interior_voxel_mesh, voxel_size(0) / 10.0, &max_min, voxel_corners);
}
/// <summary>
/// return true if I this voxel shold display in the output.
/// For the sake of this project I showed the voxels that carved (voids) as physical voxels
/// </summary>
/// <param name="vox"></param>
/// <returns></returns>
bool include_voxel_in_figure(const Voxel vox)
{
	return !vox.is_voxel_filed && vox.is_interior_voxel;
}
/// <summary>
/// triangulate matrix given the intersection of two (non self) voxels
/// </summary>
/// <param name="intersaction_vector"></param>
/// <returns></returns>
MatrixXi triangulate_matrix(std::vector<int> intersaction_vector)
{
	MatrixXi res;
	RowVector3i first_vector(intersaction_vector[0], intersaction_vector[1], intersaction_vector[2]);
	RowVector3i second_vector(intersaction_vector[2], intersaction_vector[1], intersaction_vector[3]);
	res.resize(2, 3);
	res <<
		first_vector,
		second_vector;
	return res;
}
/// <summary>
/// I don't want to show the common side of two empty voxels as well as it is not necessary to show the common side of two already displaying voxels
/// </summary>
/// <param name="first_voxel"></param>
/// <param name="second_voxel"></param>
/// <returns></returns>
bool necessary_to_append_face(const Voxel& first_voxel, const Voxel& second_voxel)
{
	if ((!include_voxel_in_figure(first_voxel) && !include_voxel_in_figure(second_voxel)) || (include_voxel_in_figure(first_voxel) && include_voxel_in_figure(second_voxel))) return true;
		return false;
}

/// <summary>
/// as input two voxels that should be neighabors, so I can make them into triangular mat and push the faces of imput faces matrix
/// </summary>
/// <param name="pushed_matrix"></param>
/// <param name="first_voxel"></param>
/// <param name="second_voxel"></param>
/// <param name="vertices_matrix"></param>
void triangulte_common_side_voxels(MatrixXi& pushed_matrix, const Voxel &first_voxel, const Voxel &second_voxel, const MatrixXd &vertices_matrix)
{
	if (necessary_to_append_face(first_voxel, second_voxel)) return;
	Voxel insideVoxel = include_voxel_in_figure(first_voxel) ? first_voxel : second_voxel;
	std::vector<int> intersaction_vector(0);//store the intersection point
	std::set_intersection(first_voxel.voxel_edge_vertices_index.begin(), first_voxel.voxel_edge_vertices_index.end(), second_voxel.voxel_edge_vertices_index.begin(), second_voxel.voxel_edge_vertices_index.end(), std::inserter(intersaction_vector, intersaction_vector.begin()));
	assert(ERROR_INTERSECT == intersaction_vector.size());
	MatrixXi triangulate_mat = triangulate_matrix(intersaction_vector); 
	MatrixXd normal_mat;
	igl::per_face_normals(vertices_matrix, triangulate_mat, normal_mat);
	assert(normal_mat.row(0).dot(normal_mat.row(1)) > 0); //if they are pointing to the same direction
	RowVector3d test_vector;
	test_vector = vertices_matrix.row(triangulate_mat(0)) - insideVoxel.center.transpose();
	bool same_direction = (normal_mat.row(0).dot(test_vector) > 0) ? true : false; //if they are pointing to the same direction
	if (same_direction) triangulate_mat << RowVector3i(intersaction_vector[2], intersaction_vector[1], intersaction_vector[0]), RowVector3i(intersaction_vector[3], intersaction_vector[1], intersaction_vector[2]);
	pushed_matrix.conservativeResize(pushed_matrix.rows() + 2, pushed_matrix.cols()); //appending smaller matrix to a bigger one
	pushed_matrix.block(pushed_matrix.rows() - 2, 0, 2, 3) = triangulate_mat;
}

void Carving_State::generate_mesh(MatrixXd& vertices_matrix, MatrixXi& faces_matrix)
{
	vertices_matrix.resize(voxel_corners.rows(), voxel_corners.cols());
	vertices_matrix = voxel_corners;
	faces_matrix.resize(0, 3);
	MatrixXi faces_matrix_of_iteration;
	Voxel voxel_of_iteration;
	int inner_loop_limit, middle_loop_limit, outer_loop_limit;
	inner_loop_limit = voxel_dimensions(2) - 1;
	middle_loop_limit = voxel_dimensions(1) - 1;
	outer_loop_limit = voxel_dimensions(0) - 1;
	for (int i = 1; i < outer_loop_limit; i++){
		for (int j = 1; j < middle_loop_limit; j++){
			for (int k = 1; k < inner_loop_limit; k++){
				faces_matrix_of_iteration.resize(0, 3); //nullifying
				voxel_of_iteration = interior_voxel_mesh(i, j, k);
				if (!include_voxel_in_figure(voxel_of_iteration)) continue; //skip this voxel since it shouldn't be included anyway
				//loops starts from '1'
				triangulte_common_side_voxels(faces_matrix_of_iteration, interior_voxel_mesh(i - 1, j, k), voxel_of_iteration, voxel_corners);
				triangulte_common_side_voxels(faces_matrix_of_iteration, interior_voxel_mesh(i + 1, j, k), voxel_of_iteration, voxel_corners);
				triangulte_common_side_voxels(faces_matrix_of_iteration, interior_voxel_mesh(i, j - 1, k), voxel_of_iteration, voxel_corners);
				triangulte_common_side_voxels(faces_matrix_of_iteration, interior_voxel_mesh(i, j + 1, k), voxel_of_iteration, voxel_corners);
				triangulte_common_side_voxels(faces_matrix_of_iteration, interior_voxel_mesh(i, j, k - 1), voxel_of_iteration, voxel_corners);
				triangulte_common_side_voxels(faces_matrix_of_iteration, interior_voxel_mesh(i, j, k + 1), voxel_of_iteration, voxel_corners);
				// Add faces to output
				int expand_faces_number = faces_matrix_of_iteration.rows();
				faces_matrix.conservativeResize(faces_matrix.rows() + expand_faces_number, faces_matrix.cols()); //expanding the matrix while keeping the information so far
				faces_matrix.block(faces_matrix.rows() - expand_faces_number, 0, expand_faces_number, faces_matrix_of_iteration.cols()) = faces_matrix_of_iteration;
	}
}	}
}
/// <summary>
/// helper function for the distance calculation of center of mass
/// </summary>
/// <param name="height_change"></param>
/// <param name="distance_between_com"></param>
/// <returns></returns>
double com_distance_helper(double height_change, double distance_between_com)
{
	double dist_squarred = std::pow(distance_between_com, 2);
	double height_squarred = std::pow(height_change, 2);
	return std::sqrt(std::max(dist_squarred, height_squarred) - std::min(dist_squarred, height_squarred));
}

/// <summary>
/// this is the function we have to minimize (using the paper energy function)
/// this function calculating the absolute value distance between com point and current iteration com. Those distances projected to the plane in order to calculate the distance
/// </summary>
/// <returns></returns>
double Carving_State::calculate_com_distance()
{
	MatrixXd center_of_mass_points_matrix;
	center_of_mass_points_matrix.resize(2, 3);
	center_of_mass_points_matrix << this->optimal_com.transpose(), this->current_com.transpose();
	VectorXd  minimum_dist_faces;
	VectorXd dist_vector;
	MatrixXd minimum_dist_points;
	MatrixXd minimum_dist_normals;
	double min, max;
	min = std::numeric_limits<double>::min();
	max = std::numeric_limits<double>::max();
	igl::signed_distance(center_of_mass_points_matrix, this->vertices_plane, this->faces_plane, igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER, min, max, dist_vector, minimum_dist_faces, minimum_dist_points, minimum_dist_normals); //calculate the distance between each com point to the ground (plane)
	double distance_com_points= (center_of_mass_points_matrix.row(0) - center_of_mass_points_matrix.row(1)).squaredNorm();
	double height_change_com_points = std::abs(dist_vector(1) - dist_vector(0));
	return (com_distance_helper(height_change_com_points, distance_com_points));
}
/// <summary>
/// halt carving when one of the following conditions are true:
/// Number of iteration exceeded the maximum fixed number that I set
/// Number of voxels that erase are exceeded the maximum number that I set 
/// point of balance and center if mass are close voxel distance from each other (paper)
/// </summary>
/// <returns></returns>
bool Carving_State::should_halt_carving()
{
	double carved_voxels_ratio;
	double dimension_product = (double)voxel_dimensions.prod();
	carved_voxels_ratio = this->number_of_voxels_to_display / dimension_product;
	//22 iterations looks enough in reasonable time from my investigation. 82% empty is also a halting condition
	return ((this->number_of_iter > 22) || (carved_voxels_ratio > 0.82) || (calculate_com_distance() < this->voxel_size.maxCoeff()));
}

/// <summary>
/// change the faces direction, given the faces matrix
/// </summary>
/// <param name="faces_matrix"></param>
void swap_faces_matrix(MatrixXi& faces_matrix)
{
	int swapper;
	for (int i = 0; i < faces_matrix.rows(); i++)
	{
		swapper = faces_matrix(i, 0);
		faces_matrix(i, 0) = faces_matrix(i, 2);
		faces_matrix(i, 2) = swapper;
	}
}

/// <summary>
/// this function helps to calculate the new center of mass (com_g) by removing or adding the additional com of two voxels.
/// </summary>
/// <param name="com_g"> global center of mass</param>
/// <param name="added_voxel"></param>
/// <param name="neig_voxel"></param>
/// <param name="vertices_matrix"></param>
/// <param name="mass"></param>
void com_step(Vector3d& com_g, const Voxel &added_voxel, const Voxel &neig_voxel, const MatrixXd &vertices_matrix, double mass)
{
	MatrixXi faces_matrix;
	faces_matrix.resize(0, 3); //setting the initial size
	triangulte_common_side_voxels(faces_matrix, added_voxel, neig_voxel, vertices_matrix);
	Vector3d current_com; //the current com is different com from the global one
	center_of_mass(vertices_matrix, faces_matrix, current_com, mass); //calculate the com with the new voxel
	if (true == include_voxel_in_figure(neig_voxel)) com_g -= current_com; //this voxel is relevant and I have to include it, i.e carve it and then distract it from global center of mass
	else com_g += current_com; //this voxel should not be carved so I have to add it to the center of mass
}

void Carving_State::carve_for_one_iteration(Vector3d& com, MatrixXd &vertices_carved_plane, MatrixXi &faces_carved_plane)
{
	this->number_of_iter = this->number_of_iter + 1;
	MatrixXd vertices_cutting_plane(3,3); //the cutting plane is vertical to the plane of the figure
	MatrixXi faces_cutting_plane(1,3);    //the cutting plane is vertical to the plane of the figure
	vertices_cutting_plane << this->optimal_com.transpose(), this->current_com.transpose(), (this->optimal_com + this->gravity_direction).transpose();
	faces_cutting_plane << 0, 1, 2;
	MatrixXd normal;
	igl::per_face_normals(vertices_cutting_plane, faces_cutting_plane, normal); //compute face normals via vertex position list, face list
	double scale_gravity;
	scale_gravity = (this->voxel_dimensions.maxCoeff() * this->voxel_size.maxCoeff() * 10) / this->gravity_direction.squaredNorm();
	double scale_vectors;
	scale_vectors = (this->voxel_dimensions.maxCoeff() * this->voxel_size.maxCoeff() * 10) / normal.row(0).squaredNorm();
	vertices_carved_plane.resize(5, 3);
	vertices_carved_plane <<this->optimal_com.transpose(), (this->optimal_com + this->gravity_direction * scale_gravity).transpose(), this->optimal_com.transpose() + normal.row(0) * scale_vectors, (this->optimal_com - this->gravity_direction * scale_gravity).transpose(), this->optimal_com.transpose() - normal.row(0) * scale_vectors; faces_carved_plane.resize(4, 3);
	faces_carved_plane << 0, 1, 2, 2, 3, 0, 0, 3, 4, 4, 1, 0;
	normal.resize(0, 0); //nullify
	igl::per_face_normals(vertices_carved_plane, faces_carved_plane, normal);//compute face normals via vertex position list, face list
	assert(normal.row(0).dot(normal.row(1)) > 0); //if they are pointing to the same direction
	assert(normal.row(1).dot(normal.row(2)) > 0); //if they are pointing to the same direction
	assert(normal.row(2).dot(normal.row(3)) > 0); //if they are pointing to the same direction
	VectorXd minimum_dist_faces;
	VectorXd dist;
	MatrixXd minimum_dist_norms;
	MatrixXd minimum_dist_points;
	igl::signed_distance(this->current_com.transpose(), vertices_carved_plane, faces_carved_plane, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist, minimum_dist_faces, minimum_dist_points, minimum_dist_norms); //computes signed distance to a mesh
	if (0.0 > dist(0)) 	swap_faces_matrix(faces_carved_plane); //should be vigilent that the plane pointing com
	std::vector<std::pair<Voxel*, double>> voxel_dist_to_plane_list; //in order not to pass by value I pass pointer to voxel, this reduce significently the time complexity
	Voxel *voxel_pointer;
	int outer_loop_limit = voxel_dimensions(0) - 1;
	int middle_loop_limit = voxel_dimensions(1) - 1;
	int inner_loop_limit = voxel_dimensions(2) - 1;
	for (int i = 1; i < outer_loop_limit; i++){
		for (int j = 1; j < middle_loop_limit; j++){
			for (int k = 1; k < inner_loop_limit; k++){
				voxel_pointer = &interior_voxel_mesh(i, j, k);
				if (!voxel_pointer->is_interior_voxel || !voxel_pointer->is_voxel_filed) continue; //this voxel should not be included
				VectorXd dist;
				VectorXd minimum_dist_faces;
				MatrixXd minimum_dist_point;
				MatrixXd minimum_dist_normal;
				igl::signed_distance(voxel_pointer->center.transpose(), vertices_carved_plane, faces_carved_plane, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist,  minimum_dist_faces, minimum_dist_point, minimum_dist_normal); //voxels distance to the cutting plane
				if (0 < dist(0)) voxel_dist_to_plane_list.push_back(std::make_pair(voxel_pointer, dist(0)));
			}
		}
	}
	struct distance_comparator
	{
		inline bool operator() (const std::pair<Voxel*, double> first_voxel_dist, const std::pair<Voxel*, double> second_voxel_dist)
		{
			return (first_voxel_dist.first < second_voxel_dist.first); //comparing the voxels using the this second pair distance
		}
	};
	std::sort(voxel_dist_to_plane_list.begin(), voxel_dist_to_plane_list.end(), distance_comparator());
	for (int i = voxel_dist_to_plane_list.size() - 1; i >= std::max(0, (int)voxel_dist_to_plane_list.size() - 1100); i--){ //carving at most 1100 voxels each iteration, start the carving operation from the most far from the voxels plane
		Voxel* voxel_of_iteration_pointer = voxel_dist_to_plane_list[i].first;
		int x_index = voxel_of_iteration_pointer->index_in_grid_s.axis_x_index;
		int y_index = voxel_of_iteration_pointer->index_in_grid_s.axis_y_index;
		int z_index = voxel_of_iteration_pointer->index_in_grid_s.axis_z_index;
		voxel_of_iteration_pointer->is_voxel_filed = false; //make this voxel as carved
		this->number_of_voxels_to_display = this->number_of_voxels_to_display + 1;
		com_step(this->current_com, *voxel_of_iteration_pointer, this->interior_voxel_mesh(x_index + 1, y_index, z_index), this->voxel_corners, this->object_mass);
		com_step(this->current_com, *voxel_of_iteration_pointer, this->interior_voxel_mesh(x_index - 1, y_index, z_index), this->voxel_corners, this->object_mass);
		com_step(this->current_com, *voxel_of_iteration_pointer, this->interior_voxel_mesh(x_index, y_index + 1, z_index), this->voxel_corners, this->object_mass);
		com_step(this->current_com, *voxel_of_iteration_pointer, this->interior_voxel_mesh(x_index, y_index - 1, z_index), this->voxel_corners, this->object_mass);
		com_step(this->current_com, *voxel_of_iteration_pointer, this->interior_voxel_mesh(x_index, y_index, z_index + 1), this->voxel_corners, this->object_mass);
		com_step(this->current_com, *voxel_of_iteration_pointer, this->interior_voxel_mesh(x_index, y_index, z_index - 1), this->voxel_corners, this->object_mass);
	}
	com = this->current_com;
}
void center_of_mass(const MatrixXd& vertices_matrix, const MatrixXi& faces_matrix, Vector3d& centers_vec, double mass)
{
	num_container nums;
	row_vecotr_container vec_container;
	centers_vec.setZero(); //setting this vector into zeros
	int faces_rows_num = faces_matrix.rows();
	for (int i = 0; i < faces_rows_num; i++) { //using Heron's
		vec_container.triangle = faces_matrix.row(i);
		vec_container.X0 = vertices_matrix.row(vec_container.triangle(0));
		vec_container.X1 = vertices_matrix.row(vec_container.triangle(1));
		vec_container.X2 = vertices_matrix.row(vec_container.triangle(2));
		vec_container.A = (vec_container.X0 - vec_container.X1);
		vec_container.B = (vec_container.X0 - vec_container.X2);
		nums.a = (vec_container.X0 - vec_container.X1).norm();
		nums.b = (vec_container.X0 - vec_container.X2).norm();
		nums.c = (vec_container.X1 - vec_container.X2).norm();
		nums.s = (nums.a + nums.b + nums.c) / 2.0;
		nums.area = std::sqrt(nums.s * (nums.s - nums.a) * (nums.s - nums.b) * (nums.s - nums.c));
		vec_container.N = vec_container.A.cross(vec_container.B).normalized();
		double helper_array[] = { vec_container.N(0), vec_container.N(1), vec_container.N(2), vec_container.X0(0), vec_container.X0(1), vec_container.X0(2), vec_container.X1(0), vec_container.X1(1), vec_container.X1(2), vec_container.X2(0), vec_container.X2(1), vec_container.X2(2) };
		centers_vec(0) += 2 * nums.area * (helper_array[0] * (helper_array[3] * helper_array[6] + helper_array[3] * helper_array[9] + helper_array[6] * helper_array[9] + std::pow(helper_array[3], 2) + std::pow(helper_array[6], 2) + std::pow(helper_array[9], 2))) / 24 + 1;
		centers_vec(1) += 2 * nums.area * (helper_array[1] * (helper_array[4] * helper_array[7] + helper_array[4] * helper_array[10] + helper_array[7] * helper_array[10] + std::pow(helper_array[4], 2) + std::pow(helper_array[7], 2) + std::pow(helper_array[10], 2))) / 24 + 1;
		centers_vec(2) += 2 * nums.area * (helper_array[2] * (helper_array[5] * helper_array[8] + helper_array[5] * helper_array[11] + helper_array[8] * helper_array[11] + std::pow(helper_array[5], 2) + std::pow(helper_array[8], 2) + std::pow(helper_array[11], 2))) / 24 + 1;
	}
	centers_vec /= mass;
}


void Carving_State::carve_whole_mesh()
{
	while(!should_halt_carving())
	{
		MatrixXd vertices_matrix_tmp;
		MatrixXi faces_matrix_tmp;
		Vector3d com_vetor_tmp;
		carve_for_one_iteration(com_vetor_tmp, vertices_matrix_tmp, faces_matrix_tmp);
	}
}
