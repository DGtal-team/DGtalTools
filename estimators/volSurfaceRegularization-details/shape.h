#pragma once

#include <DGtal/base/Common.h>
#include <DGtal/math/linalg/EigenSupport.h>

struct RoundedCubeShape
{
    typedef DGtal::Z3i::Point Point;
    typedef Eigen::Vector3d Vector;
    typedef Eigen::Matrix3d Matrix;

    RoundedCubeShape(const Matrix& transform_, const double& size_, const double& radius_) : transform(transform_), size(size_), radius(radius_)
    {
        ASSERT( radius < size );
    }

    bool operator()(const Point& point) const
    {
        const Vector vector_(point[0], point[1], point[2]);
        Vector vector = transform*vector_;
        vector = vector.cwiseAbs();
        if (vector.maxCoeff() <= size-radius) return true;
        Vector projected = vector;
        for (int kk=0; kk<projected.size(); kk++)
            if (projected[kk] > size-radius)
                projected[kk] = size-radius;
        return (vector-projected).norm() <= radius;
    }

    Matrix transform;
    double size;
    double radius;
};

struct TorusShape
{
    typedef DGtal::Z3i::Point Point;
    typedef Eigen::Vector3d Vector;
    typedef Eigen::Matrix3d Matrix;

    TorusShape(const Matrix& transform_, const double& radius_large_, const double& radius_small_) : transform(transform_), radius_large(radius_large_), radius_small(radius_small_)
    {
        ASSERT( radius_small < radius_large );
    }

    bool operator()(const Point& point) const
    {
        const Vector vector_(point[0], point[1], point[2]);
        const Vector vector = transform*vector_;
        Vector projected = vector;
        projected[1] = 0;
        if (projected.norm() == 0) return false;
        projected *= radius_large/projected.norm();
        return (vector-projected).norm() <= radius_small;
    }

    Matrix transform;
    double radius_large;
    double radius_small;
};

struct SphereShape
{
    typedef DGtal::Z3i::Point Point;
    typedef DGtal::Z3i::RealPoint RealPoint;

    SphereShape(const double& radius_, const RealPoint& center_) : radius(radius_), center(center_)
    {
    }

    bool operator()(const Point& point) const
    {
        return (point-center).norm() <= radius;
    }

    double radius;
    RealPoint center;
};

struct CapsuleShape
{
    typedef DGtal::Z3i::Point Point;
    typedef DGtal::Z3i::RealPoint RealPoint;

    CapsuleShape(const double& radius_, const double& length_, const RealPoint& direction_) : radius(radius_), length(length_), direction(direction_)
    {
        if (direction.norm() > 0) direction /= direction.norm();
    }

    bool operator()(const Point& point) const
    {
        double alpha = direction.dot(point);
        if (alpha > length/2) alpha = length/2;
        if (alpha < -length/2) alpha = -length/2;
        return (point-alpha*direction).norm() <= radius;
    }

    double radius;
    double length;
    RealPoint direction;
};

template <typename ImageType>
struct ImageShape
{
    typedef typename ImageType::Point Point;
    typedef Eigen::Vector3d Vector;
    typedef Eigen::Matrix3d Matrix;

    ImageShape(const ImageType* image_, const Point& shift_) : image(image_), shift(shift_)
    {
        ASSERT( image );
    }

    bool operator()(const Point& point_) const
    {
        const Point point = point_+shift;
        if (!image->domain().isInside(point)) return false;
        return (*image)(point) > 0;
    }

    const ImageType* image;
    Point shift;
};

struct PlaneShape
{
    typedef DGtal::Z3i::Point Point;
    typedef DGtal::Z3i::RealPoint RealPoint;

    PlaneShape(const RealPoint& normal_) : normal(normal_)
    {
    }

    bool operator()(const Point& point) const
    {
        return normal.dot(point) <= 0;
    }

    RealPoint normal;
};
