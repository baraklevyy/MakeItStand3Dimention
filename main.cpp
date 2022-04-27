#include <Eigen/Geometry>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <imgui/imgui.h>
#include <math.h>
#include "include/carving.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <iostream>
#include <vector>
#include <string>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/bounding_box.h>
#include "include/func_helper.h"
#include <igl/grid.h>
void handle_carving();
void draw_workflow_control_window();
void changing_balance_point();
bool changing_figure_orientation();
enum EXITCODE
{
	E_SUCCESS = 0,
	E_FAILED
};
void gui_vertices_update();



EXITCODE interior_mesh_clear();
EXITCODE interior_mesh_set();
int get_minimal_point_acheive(const Eigen::MatrixXd& V, bool* is_feasible);
void initGUI_current_info_s();
bool mouse_down(igl::opengl::glfw::Viewer& viewer, int x, int y);
EXITCODE refresh_mesh();
void draw_viewer_menu();
void disable_other_options(enum ACTIVE_OPERATION active_operation);
void handle_equilibrium_point();
constexpr auto ZERO_ANGLE = 0.0;
constexpr auto STRIGHT_ANGLE = 360.0;
enum ACTIVE_OPERATION
{
	POINT_SELECTION = 0,
	ANGLE_SELECTION,
	CARVING_SELECTION
};

//error code EXIT
#define Err -1
//struct to use the current identifiers to help the gui.
struct identifiers_for_View_s {
int data_interior_mesh;
int plane;
int data_mesh;
}identifiers_s;
//colors for point of balance
const Eigen::RowVector3d green(0.1,0.7,0.25);
const Eigen::RowVector3d red(0.9, 0.25, 0.333);


/*GUI variables*/
igl::opengl::glfw::imgui::ImGuiMenu menu;
igl::opengl::glfw::Viewer viewer;

//initial figures without any changes
Eigen::MatrixXi initial_F_matrix;
Eigen::MatrixXd initial_V_matrix;
//all needed data for current object
struct GUI_current_info_s
{
	Eigen::MatrixXi faces_matrix;
	//posture vector in regard to the initial figure
	Eigen::Vector3d gravity;
	Eigen::MatrixXd vertex_matrix;
	Eigen::RowVector3d center_of_plan_vec;
	int index_of_balance;
	Eigen::MatrixXi faces_plan;
	Eigen::MatrixXi interior_faces_matrix;
	Eigen::MatrixXd interior_vertex_matrix;
	Eigen::MatrixXd vertices_plan;
	//Figure position
	float pitch, yaw, roll;
	//Figure carving indication
	Carving_State *interior_mesh = NULL;
	//function selectors
	bool is_root_point_selected = true;
	bool is_figure_orientation_selected = false;
	bool is_carving = false;
} GUI_current_info_s;

//setting up the yaw, pitch and roll angles of the figure
bool changing_figure_orientation() {
	const double pi = 3.14159265358979323846;
	//we need to change the orientation; i.e change the angles of pitch, yaw and roll. We used stackoverflow for convertion euler angles to directional vector.
	double angle_convert_factor = ((pi) / (180.0));
	GUI_current_info_s.gravity(0) = -std::cos(angle_convert_factor * GUI_current_info_s.yaw) * std::sin(angle_convert_factor * GUI_current_info_s.pitch) * std::sin(angle_convert_factor * GUI_current_info_s.roll) - std::sin(angle_convert_factor * GUI_current_info_s.yaw) * std::cos(angle_convert_factor * GUI_current_info_s.roll);
	GUI_current_info_s.gravity(1) = -std::sin(angle_convert_factor * GUI_current_info_s.yaw) * std::sin(angle_convert_factor * GUI_current_info_s.pitch) * std::sin(angle_convert_factor * GUI_current_info_s.roll) + std::cos(angle_convert_factor * GUI_current_info_s.yaw) * std::cos(angle_convert_factor * GUI_current_info_s.roll);
	GUI_current_info_s.gravity(2) = std::cos(angle_convert_factor * GUI_current_info_s.pitch) * std::sin(angle_convert_factor * GUI_current_info_s.roll);
	return true;
}
void gui_vertices_update()
{
	viewer.data(identifiers_s.data_mesh).set_vertices(GUI_current_info_s.vertex_matrix);
	changing_balance_point();
}
// Update display of balance point
void disable_other_options(enum ACTIVE_OPERATION active_operation)
{
	switch (active_operation)
	{
	case POINT_SELECTION:
		GUI_current_info_s.is_carving = false;
		GUI_current_info_s.is_figure_orientation_selected = false;
		break;
	case ANGLE_SELECTION:
		GUI_current_info_s.is_carving = false;
		GUI_current_info_s.is_root_point_selected = false;
		break;
	case CARVING_SELECTION:
		GUI_current_info_s.is_root_point_selected = false;
		GUI_current_info_s.is_figure_orientation_selected = false;
		break;
	}
}
void changing_balance_point()
{
	Eigen::RowVector3d equilibrium_point = GUI_current_info_s.vertex_matrix.row(GUI_current_info_s.index_of_balance);
	viewer.data(identifiers_s.data_mesh).clear_points();
	//if we're not changing the point of balance the color of the point is green, if it on changing mode the point is red
	if (true == GUI_current_info_s.is_root_point_selected)
		viewer.data(identifiers_s.data_mesh).set_points(equilibrium_point, red);
	else
		viewer.data(identifiers_s.data_mesh).set_points(equilibrium_point, green);
}
//clearing the mesh
EXITCODE interior_mesh_clear()
{
	viewer.data(identifiers_s.data_interior_mesh).clear();
	return E_SUCCESS;
}
//settingthe mesh with matrices
EXITCODE interior_mesh_set()
{
	viewer.data(identifiers_s.data_interior_mesh).set_mesh(GUI_current_info_s.interior_vertex_matrix, GUI_current_info_s.interior_faces_matrix);
	return E_SUCCESS;
}

int get_minimal_point_acheive(const Eigen::MatrixXd& V, bool* is_feasible)
{
	double minimal_z_axis = V.col(1).minCoeff();
	int optimized_index = -500;
	for (int i = 0; i < V.rows(); i++) {
		if (V(i, 1) == minimal_z_axis) {
			optimized_index = i;
			break;
		}
	}
	if (-500 == optimized_index) *is_feasible = false;
	else
	{
		*is_feasible = true;
	}
	return optimized_index;
}
void handle_figure_angle()
{
	ImGui::Text("Please Select Figure Angle");
	if (ImGui::Checkbox("Change figure angle", &(GUI_current_info_s.is_figure_orientation_selected)))
	{
		if (true == GUI_current_info_s.is_figure_orientation_selected)
		{
			//set the angle selection option
			disable_other_options(ANGLE_SELECTION);
			//cancel the carving of the figure
			GUI_current_info_s.interior_faces_matrix.resize(0, 3);
			GUI_current_info_s.interior_vertex_matrix.resize(0, 3);
			viewer.data(identifiers_s.data_interior_mesh).clear();
			viewer.data(identifiers_s.data_mesh).show_texture = true;
			refresh_mesh();
		}
	}
	if (GUI_current_info_s.is_figure_orientation_selected)
	{
		ImGui::SliderFloat("ROLL", &(GUI_current_info_s.roll), ZERO_ANGLE, STRIGHT_ANGLE);
		if (ImGui::IsItemActive())
		{
			changing_figure_orientation();
			refresh_mesh();
		}
		ImGui::SliderFloat("PITCH", &(GUI_current_info_s.pitch), ZERO_ANGLE, STRIGHT_ANGLE);
		if (ImGui::IsItemActive())
		{
			changing_figure_orientation();
			refresh_mesh();
		}
		ImGui::SliderFloat("YAW", &(GUI_current_info_s.yaw), ZERO_ANGLE, STRIGHT_ANGLE);
		if (ImGui::IsItemActive())
		{
			changing_figure_orientation();
			refresh_mesh();
		}
	}
}
bool mouse_down(igl::opengl::glfw::Viewer& viewer, int x, int y)
{
	Eigen::RowVector3f last_mouse = Eigen::RowVector3f(
		viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0);
	if (GUI_current_info_s.is_root_point_selected)
	{
		//check for the closet mouse pointer on the figure
		int fid;
		Eigen::Vector3f bary;
		if (igl::unproject_onto_mesh(last_mouse.head(2), viewer.core().view, viewer.core().proj, viewer.core().viewport, GUI_current_info_s.vertex_matrix, GUI_current_info_s.faces_matrix, fid, bary))
		{
			long c;
			bary.maxCoeff(&c);
			GUI_current_info_s.index_of_balance = GUI_current_info_s.faces_matrix(fid, c);
			changing_balance_point();
			return true;
		}
	}
	return false;
}


void draw_viewer_menu()
{
	menu.draw_viewer_menu();
}
// move the figure obeying to the roll, pitch and yaw angles. The matrices is the conversion matrices 
EXITCODE refresh_mesh()
{
	const double pi = 3.14159265358979323846;
	const double angle_convert_factor = ((pi) / (180.0));
	const double stright_angle = 180.0;
	//setting the matrix conversion to pitch angle
	Eigen::Matrix3d pitch_matrix;
	pitch_matrix <<
		std::cos(angle_convert_factor * GUI_current_info_s.pitch), 0, -std::sin(angle_convert_factor * GUI_current_info_s.pitch),
		0, 1, 0,
		std::sin(angle_convert_factor * GUI_current_info_s.pitch), 0, std::cos(angle_convert_factor * GUI_current_info_s.pitch);
	//setting the matrix conversion to yaw angle
	Eigen::Matrix3d yaw_matrix;
	yaw_matrix <<
		std::cos(angle_convert_factor * GUI_current_info_s.yaw), -std::sin(angle_convert_factor * GUI_current_info_s.yaw), 0,
		std::sin(angle_convert_factor * GUI_current_info_s.yaw), std::cos(angle_convert_factor * GUI_current_info_s.yaw), 0,
		0, 0, 1;
	//setting the matrix conversion to roll matrix
	Eigen::Matrix3d roll_matrix;
	roll_matrix <<
		1, 0, 0,
		0, std::cos(angle_convert_factor * (double)GUI_current_info_s.roll-(stright_angle)), -std::sin(angle_convert_factor * (double)GUI_current_info_s.roll-(stright_angle)),
		0, std::sin(angle_convert_factor * (double)GUI_current_info_s.roll-(stright_angle)), std::cos(angle_convert_factor * (double)GUI_current_info_s.roll-(stright_angle));
	//always start from the initial Martix of faces and vertices.
	Eigen::Matrix3d roll_matrix_transpose = roll_matrix.transpose();
	GUI_current_info_s.vertex_matrix = initial_V_matrix * roll_matrix_transpose;
	Eigen::Matrix3d pitch_matrix_transpose = pitch_matrix.transpose();
	GUI_current_info_s.vertex_matrix = GUI_current_info_s.vertex_matrix * pitch_matrix_transpose;
	Eigen::Matrix3d yaw_matrix_transpose = yaw_matrix.transpose();
	GUI_current_info_s.vertex_matrix = GUI_current_info_s.vertex_matrix * yaw_matrix_transpose;
	Eigen::RowVector3d equilibrium_point = GUI_current_info_s.vertex_matrix.row(GUI_current_info_s.index_of_balance);
	Eigen::RowVector3d project_vector = GUI_current_info_s.center_of_plan_vec - equilibrium_point;
	Eigen::Index vertices_rows = GUI_current_info_s.vertex_matrix.rows();
	for (int i = 0; i < vertices_rows; i++) //setting the figure within the center of the plan
		GUI_current_info_s.vertex_matrix.row(i) = GUI_current_info_s.vertex_matrix.row(i) + project_vector;
	gui_vertices_update();
	return E_SUCCESS;
}


void handle_equilibrium_point()
{
	// Balance point selection
	ImGui::Text("Please Select Equilibrium Point");
	if (ImGui::Checkbox("Pick a point to balance", &(GUI_current_info_s.is_root_point_selected)))
	{
		if (true == GUI_current_info_s.is_root_point_selected) // switching to balance point selection
		{
			//disableing the other options.
			disable_other_options(POINT_SELECTION);
			changing_balance_point();
			//cancel the carving of the figure if there any
			GUI_current_info_s.interior_faces_matrix.resize(0, 3);
			GUI_current_info_s.interior_vertex_matrix.resize(0, 3);
			viewer.data(identifiers_s.data_interior_mesh).clear();
			viewer.data(identifiers_s.data_mesh).show_texture = true;
		}
	}
	if (GUI_current_info_s.is_root_point_selected)
	{
		if (ImGui::Button("Implement new point"))
		{
			//on button click: update the figure to the new point
			refresh_mesh();
		}
	}
}

void handle_carving()
{
	ImGui::Text("Carve The Figure");
	if (ImGui::Checkbox("Start Carving", &(GUI_current_info_s.is_carving)))
	{
		if (GUI_current_info_s.is_carving)
		{
			disable_other_options(CARVING_SELECTION);
			refresh_mesh();
			//here we voxel the interior mesh and make the outer voxels as transparent
			viewer.data(identifiers_s.data_mesh).show_faces = false;
			if (NULL != GUI_current_info_s.interior_mesh) delete GUI_current_info_s.interior_mesh;
			Eigen::RowVector3d equilibrium_point = GUI_current_info_s.vertex_matrix.row(GUI_current_info_s.index_of_balance);
			GUI_current_info_s.interior_mesh = new Carving_State(Eigen::Vector3d(0, -1, 0), equilibrium_point, GUI_current_info_s.vertex_matrix, GUI_current_info_s.faces_matrix, GUI_current_info_s.vertices_plan, GUI_current_info_s.faces_plan);
			GUI_current_info_s.interior_mesh->generate_mesh(GUI_current_info_s.interior_vertex_matrix, GUI_current_info_s.interior_faces_matrix);
			interior_mesh_clear();
			interior_mesh_set();
		}
	}
	if (GUI_current_info_s.is_carving)
	{
		if (ImGui::Button("Carve For One Iteration"))
		{
			if (!GUI_current_info_s.interior_mesh->should_halt_carving())
			{
				Eigen::Vector3d center_of_mass_vector;
				Eigen::MatrixXd vertices_carved_plane;
				Eigen::MatrixXi faces_carved_plane;
				GUI_current_info_s.interior_mesh->carve_for_one_iteration(center_of_mass_vector, vertices_carved_plane, faces_carved_plane);
				GUI_current_info_s.interior_mesh->generate_mesh(GUI_current_info_s.interior_vertex_matrix, GUI_current_info_s.interior_faces_matrix);
				interior_mesh_clear();
				interior_mesh_set();
			}
		}
		if (ImGui::Button("Carve To Finish"))
		{
			GUI_current_info_s.interior_mesh->carve_whole_mesh();
			GUI_current_info_s.interior_mesh->generate_mesh(GUI_current_info_s.interior_vertex_matrix, GUI_current_info_s.interior_faces_matrix);
			interior_mesh_clear();
			interior_mesh_set();
		}
		if (true == GUI_current_info_s.interior_mesh->should_halt_carving())
		{
			ImGui::Text("Process of Carving Complete");
		}
	}
}
void draw_workflow_control_window()
{
	ImGui::Begin("My Implementation Of Make It Stand", nullptr, ImGuiWindowFlags_NoSavedSettings);
	handle_equilibrium_point();
	ImGui::Separator();
	handle_figure_angle();
	ImGui::Separator();
	handle_carving();
	ImGui::End();
}
void initGUI_current_info_s()
{
	const std::string FIGURE_FILE = "../data/UrzanSoldierFiringA.stl";
	//initialize the figure file
	igl::read_triangle_mesh(FIGURE_FILE, initial_V_matrix, initial_F_matrix);
	Eigen::MatrixXd original_vertex_matrix = initial_V_matrix;
	initial_V_matrix = mesh_axises_align(original_vertex_matrix); //set the figure alined to axies
	//set the struct
	GUI_current_info_s.faces_matrix = initial_F_matrix;
	GUI_current_info_s.vertex_matrix = initial_V_matrix;

	//generate a bounding box for the figure
	Eigen::AlignedBox3d b_b;
	b_b = generate_alinged_box(initial_V_matrix);
	GUI_current_info_s.vertices_plan.resize(4, 3);
	//setting the corners (vertices) for the bounding box
	GUI_current_info_s.vertices_plan <<
		b_b.corner(b_b.BottomLeftCeil).transpose(),
		b_b.corner(b_b.BottomRightCeil).transpose(),
		b_b.corner(b_b.BottomLeftFloor).transpose(),
		b_b.corner(b_b.BottomRightFloor).transpose();
	GUI_current_info_s.faces_plan.resize(2, 3);
	GUI_current_info_s.faces_plan <<
		0, 1, 2,
		2, 1, 3;
	//getting the center of the figure plane by summing the 4 coloums and averging
	GUI_current_info_s.center_of_plan_vec = GUI_current_info_s.vertices_plan.colwise().sum() / 4.0;
	//default figure angles
	GUI_current_info_s.yaw = 0.0;
	GUI_current_info_s.pitch = 0.0;
	GUI_current_info_s.roll = 180.0;
	changing_figure_orientation();
	bool is_feasible;
	GUI_current_info_s.index_of_balance = get_minimal_point_acheive(initial_V_matrix, &is_feasible);
	if (false == is_feasible) { //point is not found
		std::cerr << "lowest mesh point is not feasible";
		exit(Err);
	}
}

int main()
{
	initGUI_current_info_s();
	viewer.plugins.push_back(&menu);
	menu.callback_draw_viewer_menu = draw_viewer_menu;
	menu.callback_draw_custom_window = draw_workflow_control_window;
	viewer.callback_mouse_down = mouse_down;
	identifiers_s.data_mesh = viewer.data().id;
	viewer.data(identifiers_s.data_mesh).set_mesh(GUI_current_info_s.vertex_matrix, GUI_current_info_s.faces_matrix);
	refresh_mesh();
	viewer.append_mesh();
	identifiers_s.data_interior_mesh = viewer.data().id;
	viewer.append_mesh();
	identifiers_s.plane = viewer.data().id;
	viewer.data(identifiers_s.plane).set_mesh(GUI_current_info_s.vertices_plan, GUI_current_info_s.faces_plan);
	viewer.selected_data_index = viewer.mesh_index(identifiers_s.data_mesh);
	viewer.core().align_camera_center(GUI_current_info_s.vertices_plan, GUI_current_info_s.faces_plan);
	viewer.launch();
	return E_SUCCESS;
}
