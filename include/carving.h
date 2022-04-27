#include <Eigen/Core>
#include "func_helper.h"
#include <set>
#include <Eigen/Geometry>
//best way to initialize 3d vector - (stackoverflow)
template <typename T>
class List3d {
public:
    List3d(size_t d1 = 0, size_t d2 = 0, size_t d3 = 0)
    {
        this->d1 = d1;
        this->d2 = d2;
        this->d3 = d3;
    }
    T& operator()(size_t i, size_t j, size_t k)
    {
        return data[i * d2 * d3 + j * d3 + k];
    }

    T const& operator()(size_t i, size_t j, size_t k) const
    {
        return data[i * d2 * d3 + j * d3 + k];
    }

    void resize(size_t d1, size_t d2, size_t d3) {
        this->d1 = d1;
        this->d2 = d2;
        this->d3 = d3;
        data.resize(d1 * d2 * d3);
    }

private:
    size_t d1, d2, d3;
    std::vector<T> data;
};


class Voxel {
public:
	bool is_voxel_filed; //if this voxel is filled or carved
	bool is_interior_voxel; //if the voxel is inside the mesh
	struct indices
	{
		int axis_x_index;
		int axis_y_index;
		int axis_z_index;
	};
	struct indices index_in_grid_s; //the indices of the voxel in the 3-Dimentional grid
	std::set<int> voxel_edge_vertices_index;
	Eigen::Vector3d center;
	double closet_dist_to_mesh; // distance to closest face on mesh
	Voxel()
	{
		is_voxel_filed = true;
		is_interior_voxel = false;
		closet_dist_to_mesh = 0;
		index_in_grid_s.axis_x_index = -1;
		index_in_grid_s.axis_y_index = -1;
		index_in_grid_s.axis_z_index = -1;
		center = Eigen::Vector3d::Zero();
	}
};

class Carving_State {
	Eigen::MatrixXi faces_plane;
	int number_of_iter;
	Eigen::Vector3d gravity_direction;
	Eigen::Vector3d current_com; 
	Eigen::Vector3i voxel_dimensions;
	double object_mass;
	int number_of_voxels_to_display;
	Eigen::Vector3d voxel_size;
	Eigen::MatrixXd voxel_corners;
	List3d<Voxel> interior_voxel_mesh;
	Eigen::Vector3d equilibrium_point;
	Eigen::MatrixXd vertices_plane; 
	Eigen::Vector3d optimal_com; 
public:
	Carving_State(const Eigen::Vector3d& gravity_direction_vector, const Eigen::Vector3d& equilibrium_point, const Eigen::MatrixXd & vertices_matrix, const Eigen::MatrixXi &faces_matrix, const Eigen::MatrixXd &vertices_plane, const Eigen::MatrixXi &faces_plane);
	double calculate_com_distance();
	void generate_mesh(Eigen::MatrixXd &vertices_matrix, Eigen::MatrixXi &faces_matrix);
	bool should_halt_carving();
	void carve_for_one_iteration(Eigen::Vector3d& com, Eigen::MatrixXd &vertices_carverd_plane, Eigen::MatrixXi &faces_carved_plane);
	void carve_whole_mesh();
};