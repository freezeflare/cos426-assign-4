// Source file for the particle system



// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
using namespace std;
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#endif




////////////////////////////////////////////////////////////
// Random Number Generator
////////////////////////////////////////////////////////////

double RandomNumber(void)
{
#ifdef _WIN32
  // Seed random number generator
  static int first = 1;
  if (first) {
    srand(GetTickCount());
    first = 0;
  }

  // Return random number
  int r1 = rand();
  double r2 = ((double) rand()) / ((double) RAND_MAX);
  return (r1 + r2) / ((double) RAND_MAX);
#else 
  // Seed random number generator
  static int first = 1;
  if (first) {
    struct timeval timevalue;
    gettimeofday(&timevalue, 0);
    srand48(timevalue.tv_usec);
    first = 0;
  }

  // Return random number
  return drand48();
#endif
}



////////////////////////////////////////////////////////////
// Generating Particles
////////////////////////////////////////////////////////////

void GenerateParticles(R3Scene *scene, double current_time, double delta_time)
{
  // Generate new particles for every source
	
	int npsources = scene -> NParticleSources();

	for (int i = 0; i < npsources; i++)
	{
		R3ParticleSource* ps = scene -> ParticleSource(i); 

		R3ShapeType type = ps -> shape -> type;
	
		if (type == R3_SPHERE_SHAPE)
		{
			R3Sphere* s = ps -> shape -> sphere;
			double r = s -> Radius();
			R3Point cen = s -> Center();
			double cen_x = cen.X();
			double cen_y = cen.Y();
			double cen_z = cen.Z();

			double z = RandomNumber() * 2 * r - r;
			double phi = RandomNumber() * 2 * M_PI;

			double d = sqrt(r * r - z * z);
			double px = cen_x + d * cos(phi);
			double py = cen_y + d * sin(phi);
			double pz = cen_z + z;

			R3Point rand_pt(px, py, pz);
			R3Vector V = rand_pt - cen;
			V *= RandomNumber();
			V *= ps -> velocity;
			
			R3Particle* p = new R3Particle();
			p -> position = rand_pt;
			p -> velocity = V;
			p -> mass = ps -> mass;
			p -> fixed = ps -> fixed;
			p -> drag = ps -> drag;
			p -> elasticity = ps -> elasticity;
			p -> lifetime = ps -> lifetime;
			p -> material = ps -> material;

			scene -> particles.push_back(p);
		}
	}

}

//Delete Particle
void DeleteParticle(R3Scene *scene, int index)
{
	int nparticles = scene -> NParticles() - 1;
	scene -> particles[index] = scene -> particles[nparticles];
	scene -> particles.pop_back();
}

//Calculate force for sink
bool CalcSinkForce(R3Scene *scene, R3Vector* f, int index, R3Point new_pos)
{
	int nsinks = scene -> NParticleSinks();
	for (int j = 0; j < nsinks; j++)
	{
		R3ParticleSink* sink = scene -> ParticleSink(j);
		double ca = sink -> constant_attenuation;
		double la = sink -> linear_attenuation;
		double qa = sink -> quadratic_attenuation;

		R3ShapeType type = sink -> shape -> type;

		if (type == R3_SPHERE_SHAPE)
		{
			R3Sphere* s = sink -> shape -> sphere;
			double r = s -> Radius();
			R3Point cen = s -> Center();

			double d = R3Distance(new_pos, cen) - r;
			if (d < 0)
			{
				DeleteParticle(scene, index);
			  return false;
			}
			else
			{
				R3Vector V = cen - new_pos;
				V.Normalize();
				V *= sink -> intensity/(ca + la * d + qa * d * d);
				*f += V;
			}
		}
	}
	return true;
}

////////////////////////////////////////////////////////////
// Updating Particles
////////////////////////////////////////////////////////////

void UpdateParticles(R3Scene *scene, double current_time, double delta_time, int integration_type)
{
  // Update position for every particle

	int nparticles = scene -> NParticles();
	R3Vector gravity = scene -> gravity;

	for (int i = 0; i < nparticles; i++)
	{
		R3Particle * cur_par = scene -> Particle(i);
		//particle lifetimes 
		if (cur_par -> lifetime < current_time && cur_par -> lifetime > 0)
		{
			DeleteParticle(scene, i);
			nparticles--;
			i--;
			continue;
		}

		R3Point old_pos = cur_par -> position;
		R3Vector velocity = cur_par -> velocity;
		R3Point new_pos = old_pos;
		new_pos += velocity * delta_time;

		R3Vector f = gravity - cur_par -> drag * velocity;

		if (!CalcSinkForce(scene, &f, i, new_pos))
		{
			nparticles--;
			i--;
			continue;
		}

		for (unsigned int k = 0; k < cur_par -> springs.size(); k++)
		{
			R3ParticleSpring* cur_spring = cur_par -> springs[k];
			R3Particle* par1;
			R3Particle* par2;
			if (cur_par == cur_spring -> particles[0])
			{
				par1 = cur_spring -> particles[0];
				par2 = cur_spring -> particles[1];
			}
			else
			{
				par2 = cur_spring -> particles[0];
				par1 = cur_spring -> particles[1];
			}
			R3Point p = par1 -> position;
			R3Point q = par2 -> position;
			R3Vector v1 = par1 -> velocity;
			R3Vector v2 = par2 -> velocity;

			R3Vector diff = q - p;
			double d = diff.Dot(diff);
			R3Vector D = diff/(d);
			double s = cur_spring -> rest_length;
			double ks = cur_spring -> ks;
			double kd = cur_spring -> kd;
			
			if (d < s)
				printf("d is %f\n", d);
			f += (ks * (d - s) + kd * (v2 - v1).Dot(D)) * D;
			//f += ks* (d - s) * D;
		}

		cur_par -> velocity += delta_time * f/ cur_par -> mass;

		cur_par -> position = new_pos;

	}
   
}



////////////////////////////////////////////////////////////
// Rendering Particles
////////////////////////////////////////////////////////////

void RenderParticles(R3Scene *scene, double current_time, double delta_time)
{
  // Draw every particle

  // REPLACE CODE HERE
  glDisable(GL_LIGHTING);
  glPointSize(5);
  glBegin(GL_POINTS);
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    glColor3d(particle->material->kd[0], particle->material->kd[1], particle->material->kd[2]);
    const R3Point& position = particle->position;
    glVertex3d(position[0], position[1], position[2]);
  }   
  glEnd();
}



