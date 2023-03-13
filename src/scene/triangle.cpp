#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL {
namespace SceneObjects {

Triangle::Triangle(const Mesh *mesh, size_t v1, size_t v2, size_t v3) {
  p1 = mesh->positions[v1];
  p2 = mesh->positions[v2];
  p3 = mesh->positions[v3];
  n1 = mesh->normals[v1];
  n2 = mesh->normals[v2];
  n3 = mesh->normals[v3];
  bbox = BBox(p1);
  bbox.expand(p2);
  bbox.expand(p3);

  bsdf = mesh->get_bsdf();
}

BBox Triangle::get_bbox() const { return bbox; }

bool Triangle::has_intersection(const Ray &r) const {
  // Part 1, Task 3: implement ray-triangle intersection
  // The difference between this function and the next function is that the next
  // function records the "intersection" while this function only tests whether
  // there is a intersection.
	//Vector3D e1(p2 - p1), e2(p3 - p1);
	//Vector3D n(cross(e1, e2));
	//if (abs(dot(n, r.d) - 0.0) <= 0.000001) return false;
	//Matrix3x3 mat(
	//	-r.d.x, e1.x, e2.x,
	//	-r.d.y, e1.y, e2.y,
	//	-r.d.z, e1.z, e2.z);
	//Vector3D tuv(mat.inv() * (r.o - p1));
	//float t(tuv.x), u(tuv.y), v(tuv.z), w(1.0 - u - v);
	//if (!((u >= 0.0) && (v >= 0.0) && (w >= 0.0))) return false;
	//if (!((t >= r.min_t) && (t <= r.max_t))) return false;
	//r.max_t = t;

  return true;

}

bool Triangle::intersect(const Ray &r, Intersection *isect) const {
  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
<<<<<<< HEAD
	//Vector3D e1(p2 - p1), e2(p3 - p1);
	//Vector3D n(cross(e1, e2));
	//if (abs(dot(n, r.d) - 0.0) <= 0.000001) return false;
	//Matrix3x3 mat(
	//	-r.d.x, e1.x, e2.x,
	//	-r.d.y, e1.y, e2.y,
	//	-r.d.z, e1.z, e2.z);
	//Vector3D tuv(mat.inv() * (r.o - p1));
	//float t(tuv.x), u(tuv.y), v(tuv.z), w(1.0 - u - v);
	//if (!((u >= 0.0) && (v >= 0.0) && (w >= 0.0))) return false;
	//if (!((t >= r.min_t) && (t <= r.max_t))) return false;
	//r.max_t = t;
	//isect->t = t;
	//isect->n = n1*w + n2*u + n3*v;
	//isect->primitive = this;
	//isect->bsdf = get_bsdf();
=======
	Vector3D e1(p2 - p1), e2(p3 - p1);
	Vector3D n(cross(e1, e2));
	if (abs(dot(n, r.d) - 0.0) <= 0.000001) return false;
	Matrix3x3 mat(
		-r.d.x, e1.x, e2.x,
		-r.d.y, e1.y, e2.y,
		-r.d.z, e1.z, e2.z);
	Vector3D tuv(mat.inv() * (r.o - p1));
	float t(tuv.x), u(tuv.y), v(tuv.z), w(1.0 - u - v);
	if (!((u >= 0.0) && (v >= 0.0) && (w >= 0.0))) return false;
	if (!((t >= r.min_t) && (t <= r.max_t))) return false;
	r.max_t = t;
	isect->t = t;
	isect->n = (n1 * w + n2 * u + n3 * v);
	isect->n.normalize();
	isect->primitive = this;
	isect->bsdf = get_bsdf();
>>>>>>> 50f0dad0cff3dde39235875bf663246bc4a7c69e

  return true;


}

void Triangle::draw(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void Triangle::drawOutline(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

} // namespace SceneObjects
} // namespace CGL
