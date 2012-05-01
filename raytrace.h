// Include file for ray tracing code

//struct for updating
struct Update_info
{
	R3Point position;
	R3Vector normal;
	double t;
};

struct R3Intersect
{
	Update_info info;
	R3Node* node;
	bool hit;
};

bool Intersect_triangle(R3Ray ray, vector<R3Point> triangle, Update_info* info);
R3Intersect Intersect_scene(R3Ray* ray, R3Scene* scene, R3Node* ignore_node);
R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth, 
  int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection, bool hard_shadow);
