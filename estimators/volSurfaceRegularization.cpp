/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file volSurfaceRegularization
 * @ingroup Tools
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * @author Pierre Gueth (\c pierre.gueth@liris.cnrs.fr )
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-smb.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 *
 *
 * @date 2011/07/07
 *
 * DGtal tool for length estimations on implicit shapes
 *
 * This file is part of the DGtal library.
 */

/**
 @page volSurfaceRegularization volSurfaceRegularization
 
 @brief Regularize a cubical complex into a smooth quadrangulated complex.

This is done by minimizing a quadratic energy function as decribed in ??. The variational
 formulation regularizes the position while aligning the regularized quads with an 
 input normal vector field. In this tool, the input normal vector field can be either 
 specified in the CSV input file, or computed using Integral Invariant (and -r option).

 @b Usage:  volSurfaceRegularization -input <volFileName> -o <regularizedcomplex.obj>

 @b Allowed @b options @b are :
 @code
 
 Positionals:
   1 TEXT REQUIRED                       input vol filename for image shape (object voxels have values > 0) or input cvs filename for surfels and normals
   1 TEXT REQUIRED                       output regularized obj

 Options:
   -h,--help                             Print this help message and exit
   -i,--image-filename TEXT REQUIRED     input vol filename for image shape (object voxels have values > 0) or input cvs filename for surfels and normals
   -o,--regularized-obj-filename TEXT REQUIRED
                                         output regularized obj
   -n,--cubical-obj-filename TEXT        output cubical obj
   -k,--shape-noise FLOAT=0              noise shape parameter
 [Option Group: Normal field estimator options]
   Options:
     -r,--normal-radius FLOAT=4            radius of normal estimator
 [Option Group: Surface approximation options]
   Options:
     -p,--regularization-position FLOAT=0.001
                                           vertex position regularization coeff
     -c,--regularization-center FLOAT=0.01 face center regularization coeff
     -a,--align FLOAT=1                    normal alignment coeff
     -f,--fairness FLOAT=0                 face fairness coeff
     -b,--barycenter FLOAT=0.1             barycenter fairness coeff
     
 @endcode

 @b Example:


 @code
 $ volSurfaceRegularization -i bunny.vol -o bunny_smooth.obj
 @endcode

 You should obtain such a result:
 @image html bunny_cubical.png "Input cubical complex."
 @image html bunny_smooth.png "Smooth quadrangulated complex."

 @see
 @ref volSurfaceRegularization.cpp

 */

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>

#include "CLI11.hpp"

#include <DGtal/images/ImageSelector.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/topology/SurfelNeighborhood.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/geometry/volumes/distance/VoronoiMap.h>
#include <random>
#include <iomanip>
#include <regex>
#include "volSurfaceRegularization-details/surface_approx.h"
#include "volSurfaceRegularization-details/surface_extract.h"
#include "volSurfaceRegularization-details/shape.h"


struct Options
{
  double noise_level {0};
  std::string image_filename;
  double normal_radius {4} ;
  double regularization_position {1e-3};
  double regularization_center {1e-2};
  double align {1.0};
  double fairness{0.0};
  double barycenter{1e-1};
  std::string regularized_obj_filename;
  std::string cubical_obj_filename;
};



  


int main(int argc, char* argv[])
{
  using DGtal::trace;
  using std::endl;
  using DGtal::PRIMAL;
  using DGtal::DUAL;
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  Options options;
  
  app.description("Regularize a cubical complex into a smooth quadrangulated complex.");
  
  app.add_option("--image-filename,-i,1",options.image_filename, "input vol filename for image shape (object voxels have values > 0) or input cvs filename for surfels and normals") -> required();
  app.add_option("--regularized-obj-filename,-o,1", options.regularized_obj_filename, "output regularized obj") -> required();
  app.add_option("--cubical-obj-filename,-n", options.cubical_obj_filename, "output cubical obj");
  app.add_option("--shape-noise,-k", options.noise_level,"noise shape parameter", true );
  auto groupNormEst = app.add_option_group("Normal field estimator options");
  groupNormEst->add_option("--normal-radius,-r", options.normal_radius, "radius of normal estimator", true);

  auto groupApprox = app.add_option_group("Surface approximation options");
  groupApprox->add_option("--regularization-position,-p", options.regularization_position, "vertex position regularization coeff", true );
  groupApprox->add_option("--regularization-center,-c",options.regularization_center, "face center regularization coeff", true);
  groupApprox->add_option("--align,-a", options.align, "normal alignment coeff", true);
  groupApprox->add_option("--fairness,-f",options.fairness,"face fairness coeff" ,true);
  groupApprox->add_option("--barycenter,-b", options.barycenter, "barycenter fairness coeff", true);

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

    const KSpace kspace;

    const Calculus calculus;
    const FlatVector original_face_normals;

    if (!ends_with(options.image_filename, ".csv"))
    {
        ASSERT( !options.image_filename.empty() );
        typedef DGtal::Z3i::Domain Domain;
        typedef DGtal::ImageSelector<Domain, unsigned char>::Type ImageUChar;
        trace.info() << "image_filename=" << options.image_filename << endl;
        ImageUChar image_uchar = DGtal::GenericReader<ImageUChar>::import(options.image_filename);
        const Domain domain = image_uchar.domain();
        trace.info() << "domain=" << domain << endl;
        const Point center = (domain.upperBound()+domain.lowerBound())/2;
        trace.info() << "center=" << center << endl;
        const ImageShape<ImageUChar> shape(&image_uchar, center);

        const_cast<KSpace&>(kspace).init(domain.lowerBound()-center-Point::diagonal(1), domain.upperBound()-center+Point::diagonal(1), true);

        std::tie(const_cast<Calculus&>(calculus), const_cast<FlatVector&>(original_face_normals)) =
            initCalculusAndNormalsWithNoise(kspace, shape, options.normal_radius, options.noise_level);
    }

    if (ends_with(options.image_filename, ".csv"))
    {
        ASSERT( !options.image_filename.empty() );
        trace.info() << "csv_filename=" << options.image_filename << endl;
        std::tie(const_cast<Calculus&>(calculus), const_cast<FlatVector&>(original_face_normals)) =
            initCalculusAndNormalsFromSurfelNormalsCSV(options.image_filename);
    }

    const FlatVector original_vertex_normals = vertexNormals(calculus, original_face_normals);

    const FlatVector original_positions;
    const FlatVector regularized_positions;
    const FlatVector original_centers;
    const FlatVector regularized_centers;
    std::tie(const_cast<FlatVector&>(original_positions),
             const_cast<FlatVector&>(regularized_positions),
             const_cast<FlatVector&>(original_centers),
             const_cast<FlatVector&>(regularized_centers)) =
             approximateSurface(calculus, original_face_normals,
             ApproxParams({options.regularization_position, options.regularization_center, options.align, options.fairness, options.barycenter}));
    ASSERT( original_positions.size() == 3*calculus.kFormLength(0, PRIMAL) );
    ASSERT( regularized_positions.size() == 3*calculus.kFormLength(0, PRIMAL) );
    ASSERT( original_centers.size() == 3*calculus.kFormLength(2, PRIMAL) );
    ASSERT( regularized_centers.size() == 3*calculus.kFormLength(2, PRIMAL) );

    {
      trace.beginBlock( "computing energies" );

      {
        double position_energy = 0;
        double align_energy    = 0;
        std::tie( position_energy, align_energy ) = approximateSurfaceEnergies(
        calculus, original_face_normals, original_positions );
        align_energy *= options.align;
        position_energy *= options.regularization_position;
        trace.info() << "original_energies=" << position_energy << " "
                     << align_energy << " " << position_energy + align_energy
                     << endl;
        }

        {
            double position_energy = 0;
            double align_energy = 0;
            std::tie(position_energy, align_energy) = approximateSurfaceEnergies(calculus, original_face_normals, regularized_positions);
            align_energy *= options.align;
            position_energy *= options.regularization_position;
            trace.info() << "regularized_energies=" << position_energy << " " << align_energy << " " << position_energy+align_energy << endl;
        }

        trace.endBlock();
    }

    {
        ASSERT( !options.regularized_obj_filename.empty() );
        trace.info() << "regularized_obj_filename=" << options.regularized_obj_filename << endl;
        exportOBJ(calculus, regularized_positions, options.regularized_obj_filename);
    }

    if (!options.cubical_obj_filename.empty())
    {
        ASSERT( !options.cubical_obj_filename.empty() );
        trace.info() << "cubical_obj_filename=" << options.cubical_obj_filename << endl;
        exportOBJ(calculus, original_positions, options.cubical_obj_filename);
    }

    return 0;
}


