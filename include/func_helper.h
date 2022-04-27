#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
/// <summary>
/// align the mesh to the axes
/// </summary>
/// <param name="vertices_matrix"></param>
/// <returns></returns>
Eigen::MatrixXd mesh_axises_align(const Eigen::MatrixXd& vertices_matrix);
/// <summary>
/// This function changing the mesh vertices from the given BB to the desired figure using the bb.
/// The calculation translation is due to the bottom left of the biunding box - Eigen library gives us this ability.
/// Maximum angle change is up to right angle
/// </summary>
/// <param name="vertices_matrix">matrix of number of vertives X 3</param>
/// <param name="desired_bounding_box">the new biunding bix the the previous matrix should fit within</param>
/// <param name="is_output_same_size_as_original">indication if we want to change the origianl size of vertices matrix</param>
/// <returns></returns>
Eigen::MatrixXd scale_and_transform_vertices(const Eigen::MatrixXd& vertices_matrix, const Eigen::AlignedBox3d &desired_bounding_box,  bool is_output_same_size_as_original);
/// <summary>
/// calculate a grid given a mesh (matrix of vertices)
/// </summary>
/// <param name="vertices_matrix"></param>
/// <param name="center_matrix"></param>
/// <param name="voxels_centers"></param>
/// <param name="grid_dimensions_in_voxel"></param>
/// <param name="voxel_size"></param>
void grid_from_mesh(
	const Eigen::MatrixXd& vertices_matrix,
	Eigen::MatrixXd& voxels_centers,       //(x + 1)(y + 1)(z + 1) X 3
	Eigen::Vector3i& grid_dimensions_in_voxel,    //the actual grid dimension in voxels 3 X 1 
	Eigen::Vector3d& voxel_size,   //voxels size 3 X 1
	Eigen::MatrixXd& center_matrix); //xyz X 3
/// <summary>
/// return box aligned to axis that calculated by vertices matrix
/// </summary>
/// <param name="vertices_matrix"></param>
/// <returns></returns>
Eigen::AlignedBox3d generate_alinged_box(const Eigen::MatrixXd& vertices_matrix);


