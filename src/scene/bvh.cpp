#include "bvh.h"
#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  BBox bbox;

  for (auto p = start; p != end; p++) {
      BBox bb = (*p)->get_bbox();
      bbox.expand(bb);
  }

  BVHNode* node = new BVHNode(bbox);

  if ((end - start) <= max_leaf_size) {

      node->start = start;
      node->end = end;

      return node;
  }

  double lowest_cost = INF_D;
  for (int i = 0; i < 3; i++) {
        for (int j = 1; j < 16; j++) {
            BBox left, right;
            int left_b(0), right_b(0);
            double plane = bbox.min[i] + ((bbox.max[i] - bbox.min[i]) / 16.0) * j;
            for (auto p = start; p != end; p++) {
                BBox bb = (*p)->get_bbox();
                if (bb.centroid()[i] < plane) {
                    left.expand(bb);
                    left_b++;
                }
                else {
                    right.expand(bb);
                    right_b++;
                }
            }
            double cost = ((double)left_b) * left.surface_area()
                + ((double)right_b)*right.surface_area();
            if (!(right_b && left_b)) continue;
            if (cost < lowest_cost) {
                axis = i; 
                bucket = j;
                lowest_cost = cost;
            }
        }
  }


  plane = bbox.min[axis] + ((bbox.max[axis] - bbox.min[axis]) / 16.0) * bucket;

  std::vector<Primitive*>::iterator bound = partition(start, end, [this](Primitive* p) {
          BBox bb = p->get_bbox();
          if (bb.centroid()[axis] < plane) {
              return true;
          }
          else {
              return false;
          }
      });
  //cout << end - start << ", " << bound - start << ", " << end - bound << endl;
  node->l = construct_bvh(start, bound, max_leaf_size);
  node->r = construct_bvh(bound, end, max_leaf_size);
  return node;

}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.

    
    double t0 = ray.min_t, t1 = ray.max_t;
    if (!(node->bb.intersect(ray, t0, t1))) return false;
    if (node->isLeaf()) {
        for (auto p = node->start; p != node->end; p++) {
            total_isects++;
            if ((*p)->has_intersection(ray))
                return true;
        }
        return false;
    }
    else {
        return (has_intersection(ray, node->l) || has_intersection(ray, node->r));
    }

}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.

    
    double t0, t1;
    if (!(node->bb.intersect(ray, t0, t1))) return false;
    if (node->isLeaf()) {
        bool hit = false;
        for (auto p = node->start; p != node->end; p++) {
            total_isects++;
            hit = (*p)->intersect(ray, i) || hit;
        }
        return hit;
    }
    else {
        bool hit1 = intersect(ray, i, node->l);
        bool hit2 = intersect(ray, i, node->r);
        return hit1 || hit2;
    }
}

} // namespace SceneObjects
} // namespace CGL
