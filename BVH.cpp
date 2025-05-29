#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod method)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(method),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    // printf( "This is construction, %d %d\n", (int)method, (int)splitMethod);
    root = recursiveBuild(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // printf( "build,%d,", (int)splitMethod );
    // exit(-1);

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim) {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x <
                       f2->getBounds().Centroid().x;
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y <
                       f2->getBounds().Centroid().y;
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z <
                       f2->getBounds().Centroid().z;
            });
            break;
        }

        std::vector<Object* >::iterator beginning, middling, ending;

        if ( splitMethod == SplitMethod::NAIVE ) {
            beginning = objects.begin();
            middling = objects.begin() + (objects.size() / 2);
            ending = objects.end();
        } else if ( splitMethod == SplitMethod::SAH ) {
            int B = 10;
            float SN = centroidBounds.SurfaceArea();
            int mincostIndex = 0;
            float minCost = std::numeric_limits<float>::infinity();
            
            for (int i = 1; i < B; i++){
                int pos = objects.size() * i / B;
                Bounds3 leftBounds, rightBounds;
                for (int k = 0; k < pos; ++k)
                    leftBounds = Union(leftBounds, objects[k]->getBounds().Centroid());
                for (int k = pos; k < objects.size(); ++k)
                    rightBounds = Union(rightBounds, objects[k]->getBounds().Centroid());
                float SA = leftBounds.SurfaceArea();
                float SB = rightBounds.SurfaceArea();
                float cost = 0.125 + (pos * SA + (objects.size() - pos) * SB) / SN;
                if (cost < minCost){
                    minCost = cost;
                    mincostIndex = i;
                }
            }
            
            beginning = objects.begin();
            middling = objects.begin() + (objects.size() * mincostIndex / B);
            ending = objects.end();
        } else {
            assert(false && "Unknown BVH split method");
        }

        auto leftshapes = std::vector<Object*>(beginning, middling);
        auto rightshapes = std::vector<Object*>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }

    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    Intersection result;
    result.happened = false;
    if ( node == nullptr )
        return result;
    auto direction_inv_neg = std::array<int, 3> {{
        int(ray.direction.x < 0),
        int(ray.direction.y < 0),
        int(ray.direction.z < 0), 
    }};

    if (!node->bounds.IntersectP(ray, ray.direction_inv, direction_inv_neg))
        return result;

    if (node->left == nullptr && node->right == nullptr) {
        result = node->object->getIntersection(ray);
        return result;
    }
    Intersection hitLeft = getIntersection(node->left, ray);
    Intersection hitRight = getIntersection(node->right, ray);
    if (!hitLeft.happened && !hitRight.happened)
        return result;
    if (!hitRight.happened || hitLeft.distance < hitRight.distance)
        return hitLeft;
    else
        return hitRight;
}