// Source file for raytracing code


// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "raytrace.h"
#include "float.h"

using namespace std;



////////////////////////////////////////////////////////////////////////
// Create image from scene
//
// This is the main ray tracing function called from raypro
// 
// "width" and "height" indicate the size of the ray traced image
//   (keep these small during debugging to speed up your code development cycle)
//
// "max_depth" indicates the maximum number of secondary reflections/transmissions to trace for any ray
//   (i.e., stop tracing a ray if it has already been reflected max_depth times -- 
//   0 means direct illumination, 1 means one bounce, etc.)
//
// "num_primary_rays_per_pixel" indicates the number of random rays to generate within 
//   each pixel during antialiasing.  This argument can be ignored if antialiasing is not implemented.
//
// "num_distributed_rays_per_intersection" indicates the number of secondary rays to generate
//   for each surface intersection if distributed ray tracing is implemented.  
//   It can be ignored otherwise.
// 
////////////////////////////////////////////////////////////////////////


//check if ray intersects a sphere, returns a bool, also updates the info
bool Intersect_ray_sphere(R3Ray ray, R3Sphere* s, R3Point p0, 
		                      Update_info* info)
{
	 p0 = ray.Start();
	 double x0 = p0.X();
	 double y0 = p0.Y();
	 double z0 = p0.Z();

	 R3Vector ray_vect = ray.Vector();
	 double x1 = ray_vect.X();
	 double y1 = ray_vect.Y();
	 double z1 = ray_vect.Z();

	 double r = s -> Radius();
	 R3Point cent = s -> Center();
	 double xc = cent.X();
	 double yc = cent.Y();
	 double zc = cent.Z();

	 double a = x1 * x1 + y1 * y1 + z1 * z1;
	 double b = 2 * ((x0 - xc) * x1 + (y0 - yc) * y1 + (z0 -zc) * z1);
	 double c = ((x0 - xc) * (x0 - xc) + (y0 - yc) * (y0 -yc)
			          + (z0 - zc) * (z0 - zc) - r * r);

	 double discri = b * b - 4 * a * c;
	 if (discri == 0)
	 {
			double time_step = (-1 * b) /(2 * a);
			info -> t = time_step;
			R3Point pos = ray.Point(time_step);
			info -> position = pos;
			R3Vector norm = cent - pos;
			norm.Normalize();
			info -> normal = norm;
			return true;
	 }
	 else
	 {
		  double time_step_1 = (-1 * b + sqrt(b * b - 4 * a * c))/(2 * a);
			double time_step_2 = (-1 * b - sqrt(b * b - 4 * a * c))/(2 * a);

			double time_step = time_step_1;
			if (time_step_1 > time_step_2)
				time_step = time_step_2;

			info -> t = time_step;
			R3Point pos = ray.Point(time_step);
			info -> position = pos;
			R3Vector norm = cent - pos;
			norm.Normalize();
			info -> normal = -norm;
			return true;
	 }

	 return false;
}

//check if the ray intersect the plane given by the 4 constants, returns the t of intersection
//if no intersection then return -1
double Intersect_plane(R3Ray ray, double a, double b, double c, double d)
{
	R3Vector ray_vec = ray.Vector();
	R3Point start = ray.Start();

	double x0 = start.X();
	double y0 = start.Y();
	double z0 = start.Z();

	double x1 = ray_vec.X();
	double y1 = ray_vec.Y();
	double z1 = ray_vec.Z();

	//parallel
	double bot = a * x1 + b * y1 + c * z1;
	if (bot == 0)
		return -1;
	double t = -1 * (a * x0 + b * y0 + c * z0 + d) /bot;

	return t;

}

//check if ray intersects the box, returns a bool and updates the info
bool Intersect_ray_box(R3Ray ray, R3Box* box, R3Point p0, 
		                      Update_info* info)
{
	double min_t = DBL_MAX;
	
	//corners of the cube
	R3Point p000 = box -> Corner(0, 0, 0);
	R3Point p100 = box -> Corner(1, 0, 0);
	R3Point p010 = box -> Corner(0, 1, 0);
	R3Point p001 = box -> Corner(0, 0, 1);
	R3Point p110 = box -> Corner(1, 1, 0);
	R3Point p101 = box -> Corner(1, 0, 1);
	R3Point p011 = box -> Corner(0, 1, 1);
	R3Point p111 = box -> Corner(1, 1, 1);

	double xmin = box -> XMin();
	double xmax = box -> XMax();
	double ymin = box -> YMin(); 
	double ymax = box -> YMax();
	double zmin = box ->ZMin();
	double zmax = box ->ZMax();

	//check if the min/maxes are reversed
	double temp;
	if (xmin > xmax)
	{
		temp = xmin;
		xmin = xmax;
		xmax = temp;
	}

	if (ymin > ymax)
	{
		temp = ymin;
		ymin = ymax;
		ymax = temp;
	}

	if (zmin > zmax)
	{
		temp = zmin;
		zmin = zmax;
		zmax = temp;
	}



	//make faces with corners
	vector<R3Plane*> faces; 
	vector<char> tables;

	faces.push_back(new R3Plane(p000, p010, p110)); //000, 010, 110, 100
	tables.push_back('z');
	faces.push_back(new R3Plane(p000, p001, p101)); //000, 001, 101, 100
	tables.push_back('y');
	faces.push_back(new R3Plane(p000, p010, p011)); //000, 010, 011, 001
	tables.push_back('x');
	faces.push_back(new R3Plane(p010, p011, p111)); //010, 011, 111, 110
	tables.push_back('y');
	faces.push_back(new R3Plane(p001, p011, p111)); //001, 011, 111, 101
	tables.push_back('z');
	faces.push_back(new R3Plane(p100, p110, p111)); //100, 110, 111, 101
	tables.push_back('x');

	int min_index = -1;

	for (unsigned int i = 0; i < faces.size(); i++)
	{
		double a = faces[i] -> A();
		double b = faces[i] -> B();
		double c = faces[i] -> C();
		double d = faces[i] -> D();

		//get the t of intersection
		double t = Intersect_plane(ray, a, b, c, d);
		//if "does not intersect"
		if (t == -1) 
		{
			//check if the ray is within the plane
			R3Point fst = ray.Point(0);
			R3Point snd = ray.Point(1);

			double sum1 = a * fst.X() + b * fst.Y() + c * fst.Z() + d;
			double sum2 = a * snd.X() + b * snd.Y() + c * snd.Z() + d;
			if (sum1 == 0 && sum2 == 0)
			{
				t = 0;
			}
			//else it does not intersect, so we skip 
			else
				continue;
		}
		R3Point intersect = ray.Point(t);
		
		double x = intersect.X();
		double y = intersect.Y();
		double z = intersect.Z();

		//check if it is inside the rectangle
		bool inside_rect = false;
		switch (tables[i])
		{
			case 'x':
				if (y > ymin && y < ymax
						&& z > zmin && z < zmax)
					inside_rect = true;
				break;
			case 'y':
				if (x > xmin && x < xmax
						&& z > zmin && z < zmax)
					inside_rect = true;
				break;
			case 'z':
				if (x > xmin && x < xmax
						&& y > ymin && y < ymax)
					inside_rect = true;
				break;
		}
		
		//update if it is
		if (t < min_t && inside_rect)
		{
			min_t = t;
			min_index = i;
		}
	}

	//if min_t is the same, then there is nothing to update
	if (min_t == DBL_MAX)
		return false;

	info -> t = min_t;
	info -> normal = faces[min_index] -> Normal();
	info -> position = ray.Point(min_t);

	return true;
}

//see if the ray intersects the triangle, returns a bool, also updates info to have the minimum t
bool Intersect_triangle(R3Ray ray, vector<R3Point> triangle, Update_info* info)
{
	R3Point p0 = ray.Start();
	R3Point p1 = triangle[0];
	R3Point p2 = triangle[1];
	R3Point p3 = triangle[2];
	
	R3Vector n = p2 - p1;
	R3Vector v2 = p3 - p1;
	n.Cross(v2);

	double a = n.X();
	double b = n.Y();
	double c = n.Z();

	double d = -(a * p1.X() + b * p1.Y() + c * p1.Z());

	double t = Intersect_plane(ray, a, b, c, d);
	if (t == -1)
		return false;
	R3Point intersect = ray.Point(t);



	for (unsigned int i = 0; i < triangle.size(); i++)
	{
		R3Vector V1 = triangle[i] - p0;
		R3Vector V2 = triangle[i + 1 != triangle.size() ? i+1 : 0] - p0;

		R3Vector N1 = V2;
		N1.Cross(V1);
		N1.Normalize();

		R3Plane* p = new R3Plane(p0, N1);
		double dist = R3SignedDistance(*p, intersect);

		if (dist < 0)
			return false;
	}

	info -> t = t;
	info -> normal = n;
	info -> position = intersect;
	return true;
}


//see if the ray intersects all the faces in the mesh, returns a bool, also updates info
bool Intersect_ray_mesh(R3Ray ray, R3Mesh* m, R3Point p0,
		                      Update_info* info)
{
	double min_t = DBL_MAX;
	int num_faces = m -> NFaces();
	for (int i = 0; i < num_faces; i++)
	{
		R3MeshFace* cur_face = m -> Face(i);
		
		//should be all triangles, check if it is
		if (cur_face -> vertices.size() != 3) 
			continue;
		//change to point array instead of R3MeshVertex
		
		vector<R3Point> triangle;
		for (int j = 0; j < 3; j++)
		{
			triangle.push_back(cur_face -> vertices[j] -> position);
		}
		Update_info new_info;	
		if (!Intersect_triangle(ray, triangle, &new_info))
			continue;

		if (new_info.t < min_t)
		{
			min_t = new_info.t;
			info -> t = new_info.t;
			info -> position = new_info.position;
			info -> normal = new_info.normal;
		}
	}

	if (min_t == DBL_MAX)
		return false;
	return true;

}

//check if ray intersects the cylinder, returns a bool and updates the info
bool Intersect_ray_cylinder(R3Ray ray, R3Cylinder* cy, R3Point p0, 
		                      Update_info* info)
{
	info -> t = DBL_MAX;
	R3Point start = ray.Start();
	R3Vector ray_vec = ray.Vector();
	R3Point cent = cy -> Center();
	double cen_x = cent.X();
	double cen_y = cent.Y();
	double cen_z = cent.Z();

	double x1 = start.X();
	double z1 = start.Z();
	double x2  = ray_vec.X();
	double z2 = ray_vec.Z();
 	
	double h = cy -> Height();
	double r = cy -> Radius();

	double a = x2 * x2 + z2 * z2;
	double b = 2 * x2 * (x1 - cen_x) + 2 * z2 * (z1 - cen_z);
	
	double c = (x1 - cen_x) * (x1 - cen_x) + (z1 - cen_z) * (z1 - cen_z) - r * r;

	double time_step_1 = -1;
	double time_step_2 = -1;
	
	time_step_1 = (-1 * b + sqrt(b * b - 4 * a * c))/(2 * a);
	time_step_2 = (-1 * b - sqrt(b * b - 4 * a * c))/(2 * a);

	R3Point p1 = ray.Point(time_step_1);
	R3Point p2 = ray.Point(time_step_2);

	//check if it is out of the bounds of the cylinder
	if (p1.Y() < cen_y - h/2.0 || p1.Y() > cen_y + h/2.0)
		time_step_1 = -1;
	if (p2.Y() < cen_y - h/2.0 || p2.Y() > cen_y + h/2.0)
		time_step_2 = -1;

	if (time_step_1 > 0)
	{
		if (time_step_2 > 0)
		{
			if (time_step_1 > time_step_2)
			{
				info -> t = time_step_2;
				info -> position = p2;
				info -> normal = p2 - cent;
				info -> normal.Normalize();
			}
			else
			{
				info -> t = time_step_1;
				info -> position = p1;
				info -> normal = p1 - cent;
				info -> normal.Normalize();
			}
		}
		else
		{
			info -> t = time_step_1;
			info -> position = p1;
			info -> normal = p1 - cent;
			info -> normal.Normalize();
		}
	}
	else if (time_step_2 > 0)
	{
		info -> t = time_step_2;
		info -> position = p2;
		info -> normal = p2 - cent;
		info -> normal.Normalize();
	}

	//check the two circle planes
	R3Vector norm1(0, 1, 0);
	R3Vector norm2(0, -1, 0);

	R3Point plane_pt1(cen_x, cen_y + h/2.0, cen_z);
	R3Point plane_pt2(cen_x, cen_y - h/2.0, cen_z);

	R3Plane* plane1 = new R3Plane(plane_pt1, norm1);
	R3Plane* plane2 = new R3Plane(plane_pt2, norm2);

	double end_cap_t1 = Intersect_plane(ray, plane1 -> A(), plane1 -> B(), 
																		plane1 -> C(), plane1 -> D());
	double end_cap_t2 = Intersect_plane(ray, plane2 -> A(), plane2 -> B(), 
																		plane2 -> C(), plane2 -> D());

	double end_cap_t = -1;
	if (end_cap_t1 > 0)
		end_cap_t = end_cap_t1;
	if (end_cap_t2 > 0 && end_cap_t2 < end_cap_t)
		end_cap_t = end_cap_t2;

	//printf("end cap t is %f\n", end_cap_t);
	//printf("current t is %f\n", info -> t);
	//printf("smaller? %d\n", end_cap_t < info-> t);

	if (end_cap_t > 0)
	{
		//printf("yes?\n");
		//check if it is inside the circle
		R3Point intersect = ray.Point(end_cap_t);
		double check = (cen_x - intersect.X()) * (cen_x - intersect.X())
				            +(cen_z - intersect.Z()) * (cen_z - intersect.Z());
		if (check < r*r && (end_cap_t < info -> t))
		{
			//printf("actually got the end cap\n");
			info -> t = end_cap_t;
			info -> position = intersect;
			if (end_cap_t == end_cap_t1)
				info -> normal = norm1;
			else 
				info -> normal = norm2;
			return true;
		}
	}

	if (info -> t < DBL_MAX)
		return true;
	return false;
}

//check if ray intersects the cone, returns a bool and updates the info
bool Intersect_ray_cone(R3Ray ray, R3Cone* cone, R3Point p0, 
		                      Update_info* info)
{
	info -> t = DBL_MAX;
	R3Point start = ray.Start();
	R3Vector ray_vec = ray.Vector();
	R3Point cent = cone -> Center();
	double cen_x = cent.X();
	double cen_y = cent.Y();
	double cen_z = cent.Z();

	double x1 = start.X();
	double z1 = start.Z();
	double y1 = start.Y();
	double x2 = ray_vec.X();
	double z2 = ray_vec.Z();
	double y2 = ray_vec.Y();
	
	double h = cone -> Height();
	double r = cone -> Radius();
	double k = r/h;
	double y_max = cen_y + h/2.0;

	double a = x2 * x2 + z2 * z2 - k * k * y2 * y2;
	double b = 2 *(x2 * (x1 - cen_x) + z2 * (z1 - cen_z) - k * k * y2 * (y1 - y_max));
	
	double c = (x1 - cen_x) * (x1 - cen_x) + (z1 - cen_z) 
						* (z1 - cen_z) - k * k * (y1 - y_max) * (y1 - y_max);

	double time_step_1 = -1;
	double time_step_2 = -1;
	
	time_step_1 = (-1 * b + sqrt(b * b - 4 * a * c))/(2 * a);
	time_step_2 = (-1 * b - sqrt(b * b - 4 * a * c))/(2 * a);

	R3Point p1 = ray.Point(time_step_1);
	R3Point p2 = ray.Point(time_step_2);

	//check if it is out of the bounds of the cone
	if (p1.Y() < cen_y - h/2.0 || p1.Y() > cen_y + h/2.0)
		time_step_1 = -1;
	if (p2.Y() < cen_y - h/2.0 || p2.Y() > cen_y + h/2.0)
		time_step_2 = -1;

	if (time_step_1 > 0)
	{
		if (time_step_2 > 0)
		{
			if (time_step_1 > time_step_2)
			{
				info -> t = time_step_2;
				info -> position = p2;
				info -> normal = p2 - cent;
				info -> normal.Normalize();
			}
			else
			{
				info -> t = time_step_1;
				info -> position = p1;
				info -> normal = p1 - cent;
				info -> normal.Normalize();
			}
		}
		else
		{
			info -> t = time_step_1;
			info -> position = p1;
			info -> normal = p1 - cent;
			info -> normal.Normalize();
		}
	}
	else if (time_step_2 > 0)
	{
		info -> t = time_step_2;
		info -> position = p2;
		info -> normal = p2 - cent;
		info -> normal.Normalize();
	}

	//check the bot circle planes
	R3Vector norm1(0, -1, 0);

	R3Point plane_pt1(cen_x, cen_y - h/2.0, cen_z);

	R3Plane* plane1 = new R3Plane(plane_pt1, norm1);

	double end_cap_t = Intersect_plane(ray, plane1 -> A(), plane1 -> B(), 
																		plane1 -> C(), plane1 -> D());

	if (end_cap_t > 0)
	{
		R3Point intersect = ray.Point(end_cap_t);
		double check = (cen_x - intersect.X()) * (cen_x - intersect.X())
				            +(cen_z - intersect.Z()) * (cen_z - intersect.Z());
		if (check < r*r && (end_cap_t < info -> t))
		{
			info -> t = end_cap_t;
			info -> position = intersect;
			info -> normal = norm1;
			return true;
		}
	}

	if (info -> t < DBL_MAX)
		return true;
	return false;
}


//the overall fucntion that calls different functions for different shape types
bool DoIntersect(R3Ray ray, R3Node* node, Update_info* cur_info, R3Point p0)
{
	R3Shape* shape = node -> shape;

	if (shape == NULL)
		return false;
	if (shape -> type == R3_BOX_SHAPE)
	{
		return Intersect_ray_box(ray, shape -> box, p0, cur_info);
	}
	else if (shape -> type == R3_SPHERE_SHAPE)
	{
		return Intersect_ray_sphere(ray, shape -> sphere, p0, cur_info);
	}
	else if (shape -> type == R3_MESH_SHAPE)
	{
		return Intersect_ray_mesh(ray, shape -> mesh, p0, cur_info);
	}
	else if (shape -> type == R3_CYLINDER_SHAPE)
		return Intersect_ray_cylinder(ray, shape -> cylinder, p0, cur_info);
	else if (shape -> type == R3_CONE_SHAPE)
		return Intersect_ray_cone(ray, shape -> cone, p0, cur_info);

	return false;
}

//recursive function to go through a scene
R3Intersect Traverse_scene(R3Ray ray, R3Node* node, R3Intersect cur_best, R3Point p0, R3Node* ignore_node)
{
	R3Ray orig_ray = ray;
	R3Vector orig_vector = orig_ray.Vector();
	R3Point orig_start = orig_ray.Start();
	double orig_x = orig_start.X();
	double orig_y = orig_start.Y();
	double orig_z = orig_start.Z();

	vector<R3Node*> cur_child = node -> children;
	R3Matrix trans = node -> transformation;
	R3Matrix trans_inv = trans.Inverse();
	ray.Transform(trans_inv);

	//end case
	if (cur_child.size() == 0)
	{
		return cur_best;
	}

	for(unsigned int i = 0; i < cur_child.size(); i++)
	{
     if (cur_child[i] == ignore_node)
            continue;
	  
		 //bbox code, but somehow does not improve performance
		 /*Update_info bbox_info;

		if (Intersect_ray_box(ray, &(cur_child[i] -> bbox), p0, &bbox_info))
		{

			if (bbox_info.t > cur_best.info.t)
				continue;
		}*/
		

			
		Update_info cur_info;



		if (DoIntersect(ray, cur_child[i], &cur_info, p0))
		{
			double new_t = 0;
			cur_info.position.Transform(trans);
			cur_info.normal.Transform(trans);
			//calculate new t
			if (orig_x != 0)
				new_t = (cur_info.position.X() - orig_x)/orig_vector.X();
			else if (orig_y != 0)
				new_t = (cur_info.position.Y() - orig_y)/orig_vector.Y();
			else 
				new_t = (cur_info.position.Z() - orig_z)/orig_vector.Z();

			//update the minimum info
			if (new_t < cur_best.info.t && new_t > 0)
			{
				cur_best.info.t = new_t; 
				cur_best.info.position = cur_info.position;
				cur_best.info.normal = cur_info.normal;
				cur_best.info.normal.Normalize();
				

				cur_best.node = cur_child[i];
			}
		}
		cur_best = Traverse_scene(ray, cur_child[i], cur_best, p0, ignore_node);
	}
	return cur_best;
}

//function that calls the recursive scene traversal function
R3Intersect Intersect_scene(R3Ray* ray, R3Scene* scene, R3Node* ignore_node)
{
	R3Node* root = scene -> Root();

	R3Intersect max;
	max.info.t = DBL_MAX;
	R3Point p0 = scene -> camera.eye;


	R3Intersect new_min_info = Traverse_scene(*ray, root, max, p0, ignore_node);
	if (new_min_info.info.t == DBL_MAX)
		new_min_info.hit = false;
	else
		new_min_info.hit = true;


	return new_min_info;
}

R3Vector Random_vector(R3Light* light, R3Point intersect)
{
	double x = 1.0 *rand()/RAND_MAX;
	double y = 1.0 *rand()/RAND_MAX;
	double z = 1.0 *rand()/RAND_MAX;

	R3Point rand_pt(x,y,z);
	R3Vector rand_vec = light -> position - rand_pt;
	
	rand_vec.Cross(light -> direction);
	rand_vec.Normalize();

	double r = 1.0 *rand()/RAND_MAX;
	rand_vec *= r;
	R3Point copy_pos(light -> position);
	copy_pos.Translate(rand_vec);
	R3Vector ret_vec = copy_pos - intersect;
	return ret_vec;
}

//returns the Phong illumination
R3Rgb Phong(R3Scene* scene, R3Ray* ray, R3Intersect* info, double u, double v,
		        int max_depth, int cur_depth, bool hard_shadow)
{
	//printf("1\n");
	R3Node* node = info -> node;
	R2Image* texture = node -> material -> texture;
	bool isTexture = false;
	R2Pixel pix;
	if (texture != NULL)
	{
 		pix = texture -> Pixel(round(u), round(v));
		isTexture = true;
	}

	R3Rgb kd = node -> material -> kd;
	R3Rgb ka = node -> material -> ka;
	R3Rgb ks = node -> material -> ks;
	R3Rgb kt = node -> material -> kt;
	double shininess = node -> material -> shininess;
	R3Vector N = info -> info.normal;
	N.Normalize();
	R3Point intersect = info -> info.position;

	R3Point start = ray -> Start();
	R3Vector V = start - intersect;
	V.Normalize();
	R3Rgb radiance;

	radiance = ka * scene -> ambient;
	R3Rgb emission = info -> node -> material -> emission;
	//R3Rgb reflection(0,0,0);


	int NLights = scene -> NLights();
	for (int i = 0; i < NLights; i++)
	{
		R3Light* cur_light = scene -> Light(i);
		R3Rgb I0;
		R3Point light_pos = cur_light -> position;
		R3Vector D = cur_light -> direction;
		D.Normalize();
		R3Vector L;
		int loop_count = 0;
		if (cur_light -> type == R3_AREA_LIGHT)
			loop_count = 16;
		else
			loop_count = 1;
		int max_count = loop_count;
		//softshadow loop
		while (loop_count > 0)
		{
			if (cur_light -> type == R3_DIRECTIONAL_LIGHT)
				L =  -cur_light -> direction;
			else if (cur_light -> type == R3_POINT_LIGHT)
			{
				L = cur_light -> position - intersect;
			}
			else if (cur_light -> type == R3_SPOT_LIGHT)
			{
				L =  cur_light -> position - intersect;
			}
			else if (cur_light -> type == R3_AREA_LIGHT)
			{
				L = Random_vector(cur_light, intersect);
			}
			L.Normalize();
			R3Vector R = 2.0 * (L.Dot(N)) * (N) - L;

			//specular reflection
			bool need_reflect = false;
			R3Rgb ref_col;
			if (cur_depth < max_depth)
			{
				R3Ray* reflect_ray = new R3Ray(intersect, R); 
				R3Intersect reflect_info = Intersect_scene(reflect_ray, scene, node);
				if (reflect_info.hit == true)
				{
					ref_col = Phong(scene, reflect_ray, 
							&reflect_info, u, v, max_depth, cur_depth + 1, hard_shadow);
					need_reflect = true;
				}
			}

			R3Vector copy_L(L);
			copy_L.Flip();
			//transmission
			R3Rgb trans_col;
			if (cur_depth < max_depth)
			{
				R3Ray* trans_ray = new R3Ray(intersect, copy_L);
				R3Intersect trans_info = Intersect_scene(trans_ray, scene, node);
				if (trans_info.hit == true)
				{
					trans_col = Phong(scene, trans_ray, 
							&trans_info, u, v, max_depth, cur_depth + 1, hard_shadow);
				}
			}


			if (need_reflect)
				I0 += ks * ref_col;

			I0 += kt * trans_col;



			//cast hard shadow
			if (hard_shadow)
			{
				R3Ray* shadow = new R3Ray(intersect, light_pos);
				R3Intersect cur_best;
				cur_best.info.t = DBL_MAX;
				R3Intersect ret_info;
				ret_info = Traverse_scene(*shadow, scene -> root, cur_best, 
						shadow->Start(), node);


				if(ret_info.info.t < DBL_MAX && ret_info.info.t > 0)
				{
					//printf("t is %f\n", ret_info.info.t);
					double new_t = shadow->T(light_pos);
					if (ret_info.info.t < new_t)
					{
						loop_count--;  
						continue;
					}
				}
			}

			if (N.Dot(L) > 0)
			{
				if (isTexture)
					I0 += pix * kd  * (N.Dot(L));
				else
					I0 += kd  * (N.Dot(L));

			}
			if (V.Dot(R) > 0)
			{
				I0 += ks * pow(V.Dot(R), shininess);
			}


			double sc = cur_light -> angle_cutoff;
			double sd = cur_light -> angle_attenuation;

			double ca = cur_light -> constant_attenuation;
			double la = cur_light -> linear_attenuation;
			double qa = cur_light -> quadratic_attenuation;
			double d = R3Distance(intersect, light_pos);
			if (cur_light -> type == R3_POINT_LIGHT)
			{
				I0 *= 1.0/(ca + la * d + qa * d * d);
			}
			else if (cur_light -> type == R3_SPOT_LIGHT)
			{
				double theta = acos((-L).Dot(D));
				if (theta > sc)
					I0 *= 0;
				else
					I0 *= pow(cos(theta), sd)/(ca + la * d + qa * d * d);
			}
			else if (cur_light -> type == R3_AREA_LIGHT)
			{
				L.Flip();
				double dot = D.Dot(L);
				if (dot > 0)
					I0 *= dot/(ca + la * d + qa * d * d);
			}

			I0 *= cur_light -> color;
			loop_count--; 
		} 
		I0 /= max_count; 

		
		radiance += I0;
	}
	
	radiance += emission;
	return radiance;
}

R3Rgb ComputeRadiance(R3Scene* scene, R3Ray* ray, R3Intersect* info, double u, double v,
		                  int max_depth, bool hard_shadow)
{
	//R3Ray specular_ray = Specular_ray(ray, info);
	//R3Ray refractive_ray = Refractive_ray(ray, info);

	R3Rgb radiance = Phong(scene, ray, info, u, v, max_depth, 0, hard_shadow);

	return radiance;
}


R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth,
  int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection,
	bool hard_shadow)
{
  // Allocate  image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }

	R3Point p0 = scene -> camera.eye;
	R3Vector right = scene -> camera.right;
	R3Vector towards = scene -> camera.towards;
	R3Vector up = scene -> camera.up;
	double xfov = scene -> camera.xfov;
	double yfov = scene -> camera.yfov;
	
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			double x_weight = (j - width/2.0)/(width/2.0);
			double y_weight = (i - height/2.0)/(height/2.0);

			//printf("x: %d, y: %d\n", j, i);

			R3Vector vec = towards + x_weight * tan(xfov) * right + y_weight * tan(yfov) * up;
			R3Ray* ray = new R3Ray(p0, vec);


			R3Intersect info = Intersect_scene(ray, scene, NULL);
			if (info.hit)
			{
				//for texture mapping
				double u = asin(info.info.normal.X()/M_PI) + 0.5;
				double v = asin(info.info.normal.Y()/M_PI) + 0.5;

				u *= width;
				v *= height;

				R3Rgb radiance = ComputeRadiance(scene, ray, &info, u, v, max_depth, hard_shadow);
				//printf("%f %f %f", radiance[0], radiance[1], radiance[2]);
				R2Pixel* new_pix = new R2Pixel(radiance);
				image -> Pixel(j,i) = *new_pix; 
			}
			else
				image -> Pixel(j,i) = scene -> background;
		}
	}
  // Return image
  return image;
}
