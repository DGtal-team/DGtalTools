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
File options:
  -i [ --image-filename ] arg           input vol filename for image shape (object voxels
                                        have values > 0) or input cvs filename for surfels and
                                        normals
  -o [ --regularized-obj-filename ] arg output regularized obj
  -n [ --cubical-obj-filename ] arg     output cubical obj
  -k [ --shape-noise ] arg (=0)         noise shape parameter

Normal field estimator options:
  -r [ --normal-radius ] arg (=4) radius of normal estimator

Surface approximation options:
  -p [ --regularization-position ] arg (=1e-3)
                                        vertex position regularization coeff
  -c [ --regularization-center ] arg (=1e-2)
                                        face center regularization coeff
  -a [ --align ] arg (=1)               normal alignment coeff
  -f [ --fairness ] arg (=0)            face fairness coeff
  -b [ --barycenter ] arg (=1e-1)       barycenter fairness coeff
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
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/errors.hpp>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/io/readers/GenericReader.h>
#include <random>
#include <iomanip>
#include <regex>
#include "volSurfaceRegularization-details/surface_approx.h"
#include "volSurfaceRegularization-details/surface_extract.h"
#include "volSurfaceRegularization-details/shape.h"


struct Options
{
    double noise_level;
    std::string image_filename;
    double normal_radius;
    double regularization_position;
    double regularization_center;
    double align;
    double fairness;
    double barycenter;
    std::string regularized_obj_filename;
    std::string cubical_obj_filename;
};

Options
parse_options(int argc, char* argv[])
{
    namespace po = boost::program_options;

    using DGtal::trace;
    using std::endl;

    Options options;
    po::options_description po_shape("File options");
    po_shape.add_options()
        ("image-filename,i", po::value<std::string>(&options.image_filename)->default_value(""), "input vol filename for image shape (object voxels have values > 0) or input cvs filename for surfels and normals")
        ("regularized-obj-filename,o", po::value<std::string>(&options.regularized_obj_filename)->default_value(""), "output regularized obj")
        ("cubical-obj-filename,n", po::value<std::string>(&options.cubical_obj_filename)->default_value(""), "output cubical obj")
        ("shape-noise,k", po::value<double>(&options.noise_level)->default_value(0), "noise shape parameter")
        ;

    po::options_description po_normal("Normal field estimator options");
    po_normal.add_options()
        ("normal-radius,r", po::value<double>(&options.normal_radius)->default_value(4), "radius of normal estimator")
        ;

    po::options_description po_approx("Surface approximation options");
    po_approx.add_options()
        ("regularization-position,p", po::value<double>(&options.regularization_position)->default_value(1e-3, "1e-3"), "vertex position regularization coeff")
        ("regularization-center,c", po::value<double>(&options.regularization_center)->default_value(1e-2, "1e-2"), "face center regularization coeff")
        ("align,a", po::value<double>(&options.align)->default_value(1), "normal alignment coeff")
        ("fairness,f", po::value<double>(&options.fairness)->default_value(0), "face fairness coeff")
        ("barycenter,b", po::value<double>(&options.barycenter)->default_value(1e-1, "1e-1"), "barycenter fairness coeff")
        ;

    po::options_description po_options("surfaceApprox [options]");
    po_options.add(po_shape).add(po_normal).add(po_approx).add_options()
        ("help,h", "display this message")
        ;

    po::positional_options_description positional;
    positional.add("image-filename",1);
    positional.add("regularized-obj-filename",1);

    try
    {
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(po_options).positional(positional).run(), vm);
        po::notify(vm);

        if (vm.count("help"))
        {
            trace.info() << po_options;
            std::exit(0);
        }

        if (options.image_filename.empty()) throw po::validation_error(po::validation_error::invalid_option_value);
        if (options.regularized_obj_filename.empty()) throw po::validation_error(po::validation_error::invalid_option_value);
    }
    catch (std::exception& ex)
    {
        trace.error() << ex.what() << endl;
        trace.info() << po_options;
        std::exit(1);
    }

    return options;
}


int main(int argc, char* argv[])
{
    using DGtal::trace;
    using std::endl;
    using DGtal::PRIMAL;
    using DGtal::DUAL;

    const Options options = parse_options(argc, argv);

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


