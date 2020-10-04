#include "surface_approx.h"
#include <regex>
#include <DGtal/dec/DiscreteExteriorCalculusFactory.h>

bool ends_with(const std::string& value, const std::string& ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

template <typename Value>
Value
format_to(const std::string& str)
{
    std::stringstream ss;
    ss << str;
    Value foo;
    ss >> foo;
    return foo;
}

std::tuple<Calculus, FlatVector>
initCalculusAndNormalsFromSurfelNormalsCSV(const std::string& filename)
{
    using DGtal::trace;
    using std::endl;
    using DGtal::PRIMAL;

    trace.beginBlock("loading surfels and normals from csv file");

    trace.info() << "filename=" << filename << endl;

    typedef DGtal::Z3i::RealVector RealVector;
    typedef std::map<SCell, RealVector> Normals;
    typedef std::vector<SCell> SCells;

    const KSpace kspace;

    static const std::string float_re = "([-+]?([[:digit:]]*\\.?[[:digit:]]+|[[:digit:]]+\\.?[[:digit:]]*)([eE][-+]?[[:digit:]]+)?)";
    static const std::regex line_re("^ *(-?[[:digit:]]+) +(-?[[:digit:]]+) +(-?[[:digit:]]+) +([0-1]) +"+float_re+" +"+float_re+" +"+float_re+" *$", std::regex::extended);

    Normals normals;
    SCells cells;

    try
    {
        std::ifstream handle(filename);

        while (true)
        {
            std::string line;
            std::getline(handle, line);
            if (!handle.good()) break;
            std::smatch what;
            if (std::regex_match(line, what, line_re) <= 0)
                throw std::ios_base::failure("bad csv format");

            ASSERT( what.size() == 14 );

            const std::string xx_string(what[1].first, what[1].second);
            const std::string yy_string(what[2].first, what[2].second);
            const std::string zz_string(what[3].first, what[3].second);
            const std::string sign_string(what[4].first, what[4].second);
            const std::string nx_string(what[5].first, what[5].second);
            const std::string ny_string(what[8].first, what[8].second);
            const std::string nz_string(what[11].first, what[11].second);

            Point point;
            point[0] = format_to<Point::Coordinate>(xx_string);
            point[1] = format_to<Point::Coordinate>(yy_string);
            point[2] = format_to<Point::Coordinate>(zz_string);

            const bool positive = format_to<bool>(sign_string);

            const SCell cell = kspace.sCell(point, positive);

            cells.push_back(cell);

            RealPoint normal;
            normal[0] = format_to<RealPoint::Coordinate>(nx_string);
            normal[1] = format_to<RealPoint::Coordinate>(ny_string);
            normal[2] = format_to<RealPoint::Coordinate>(nz_string);

            normals.insert(std::make_pair(cell, normal));
        }
    }
    catch (std::ios_base::failure &fail)
    {
        trace.error() << fail.what() << endl;
        trace.endBlock();
        exit(42);
    }

    trace.info() << "cells_size=" << cells.size() << endl;
    trace.info() << "normals_size=" << normals.size() << endl;

    ASSERT( normals.size() == cells.size() );

    Calculus calculus = DGtal::DiscreteExteriorCalculusFactory<Backend, int32_t>::createFromNSCells<2>(cells.begin(), cells.end(), true);

    const auto buildFlatVector = [&kspace, &calculus](const Normals& real_vectors)
    {
        ASSERT( real_vectors.size() == calculus.kFormLength(2, PRIMAL) );
        const int nsurfels = (const int)real_vectors.size();
        FlatVector vectors(3*real_vectors.size());
        for (const std::pair<SCell, RealVector>& cell_pair : real_vectors)
        {
            const int index = (const int)calculus.getCellIndex(kspace.unsigns(cell_pair.first));
            for (int dim=0; dim<3; dim++)
                vectors[index+nsurfels*dim] = cell_pair.second[dim];
        }
        return vectors;
    };

    trace.endBlock();

    return std::make_tuple(calculus, buildFlatVector(normals));
}

void
exportOBJ(const Calculus& calculus, const FlatVector& positions, const std::string& filename)
{
    using DGtal::PRIMAL;
    using std::endl;

    ASSERT( !filename.empty() );

    const KSpace& kspace = calculus.myKSpace;
    const int nvertices = (const int)calculus.kFormLength(0, PRIMAL);

    ASSERT( positions.size() == 3*nvertices );
    typedef Calculus::Index Index;

    std::ofstream handle(filename);

    for (int index=0; index<nvertices; index++)
    {
        handle << "v";
        for (int dim=0; dim<3; dim++)
            handle << " " << positions[index+dim*nvertices];
        handle << endl;
    }

    for (const Calculus::ConstIterator::value_type& pair : calculus)
    {
        const Cell& face = pair.first;
        if (kspace.uDim(face) != 2) continue;

        std::vector<Cell> vertices;
        {
            std::set<Cell> vertices_set;
            for (const Cell& edge : kspace.uLowerIncident(face))
                for (const Cell& vertex : kspace.uLowerIncident(edge))
                    if (calculus.containsCell(vertex))
                        vertices_set.insert(vertex);
            ASSERT( vertices_set.size() == 4 );
            std::copy(vertices_set.begin(), vertices_set.end(), std::back_inserter(vertices));
            std::swap(vertices[2], vertices[3]);
        }

        if (pair.second.flipped)
        {
            std::swap(vertices[0], vertices[1]);
            std::swap(vertices[2], vertices[3]);
        }

        handle << "f";
        for (const Cell& cell : vertices)
        {
            const Index index = calculus.getCellIndex(cell);
            ASSERT( index >= 0 && index < nvertices );
            handle << " " << (index+1);
        }
        handle << endl;
    }
}

FlatVector
vertexNormals(const Calculus& calculus, const FlatVector& face_normals)
{
    using DGtal::PRIMAL;

    const KSpace& kspace = calculus.myKSpace;
    const int nfaces = (const int)calculus.kFormLength(2, PRIMAL);
    const int nvertices = (const int)calculus.kFormLength(0, PRIMAL);

    ASSERT( face_normals.size() == 3*nfaces );

    typedef Calculus::Index Index;

    FlatVector vertex_normals = FlatVector::Zero(3*nvertices);
    {
        Index face_index = 0;
        for (const SCell& signed_face : calculus.getIndexedSCells<2, PRIMAL>())
        {
            const Cell& face = kspace.unsigns(signed_face);

            std::set<Cell> face_vertices;
            for (const Cell& edge : kspace.uLowerIncident(face))
                for (const Cell& vertex : kspace.uLowerIncident(edge))
                    if (calculus.containsCell(vertex))
                        face_vertices.insert(vertex);
            ASSERT( face_vertices.size() == 4 );

            for (const Cell& vertex : face_vertices)
            {
                const Index vertex_index = calculus.getCellIndex(vertex);
                for (int dim=0; dim<3; dim++)
                    vertex_normals[vertex_index+dim*nvertices] += face_normals[face_index+dim*nfaces];
            }

            face_index++;
        }
    }

    for (Index vertex_index=0; vertex_index<nvertices; vertex_index++)
    {
        double norm = 0;
        for (int dim=0; dim<3; dim++)
        {
            const double foo = vertex_normals[vertex_index+dim*nvertices];
            norm += foo*foo;
        }
        norm = sqrt(norm);
        if (norm == 0) continue;
        for (int dim=0; dim<3; dim++)
            vertex_normals[vertex_index+dim*nvertices] /= norm;
    }

    return vertex_normals;
}

bool
checkOperatorSymmetry(const OperatorMatrix& matrix, const double tol)
{
    typedef std::vector<Triplet> Triplets;

    Triplets matrix_transpose_triplets;
    for (int kk=0; kk<matrix.outerSize(); kk++)
        for(OperatorMatrix::InnerIterator ii(matrix, kk); ii; ++ii)
            matrix_transpose_triplets.push_back(Triplet(ii.col(), ii.row(), ii.value()));

    OperatorMatrix matrix_transpose(matrix.cols(), matrix.rows());
    matrix_transpose.setFromTriplets(matrix_transpose_triplets.begin(), matrix_transpose_triplets.end());
    const OperatorMatrix foo = matrix - matrix_transpose;

    for (int kk=0; kk<foo.outerSize(); kk++)
        for(OperatorMatrix::InnerIterator ii(foo, kk); ii; ++ii)
            if (std::abs(ii.value()) > tol)
                return false;

    return true;
}

std::tuple<double, double>
approximateSurfaceEnergies(const Calculus& calculus, const FlatVector& normals, const FlatVector& positions)
{
    using DGtal::PRIMAL;

    typedef Calculus::Index Index;

    const KSpace& kspace = calculus.myKSpace;
    const CellEmbedder embedder(kspace);
    const int nfaces = (const int)calculus.kFormLength(2, PRIMAL);
    const int nvertices = (const int)calculus.kFormLength(0, PRIMAL);

    ASSERT( normals.size() == 3*nfaces );
    ASSERT( positions.size() == 3*nvertices );

    double position_energy = 0;
    for (const Calculus::ConstIterator::value_type& pair : calculus)
    {
        const Cell& vertex = pair.first;
        if (kspace.uDim(vertex) != 0) continue;
        const Index vertex_index = calculus.getCellIndex(vertex);

        ASSERT( pair.second.primal_size == 1 );
        const double area = pair.second.dual_size;

        const RealPoint original_position = embedder(vertex);
        RealPoint position;
        for (int dim=0; dim<3; dim++)
            position[dim] = positions[vertex_index+dim*nvertices];

        for (int dim=0; dim<3; dim++)
            position_energy += area*pow(position[dim]-original_position[dim], 2);
    }

    double align_energy = 0;
    for (const Calculus::ConstIterator::value_type& pair : calculus)
    {
        const Cell& face = pair.first;
        if (kspace.uDim(face) != 2) continue;

        ASSERT( pair.second.dual_size == 1 );
        const double area = pair.second.primal_size;
        ASSERT( area == 1 ); // true on cubical complexes

        RealPoint face_normal;
        for (int dim=0; dim<3; dim++)
            face_normal[dim] = normals[pair.second.index+dim*nfaces];

        std::vector<Cell> vertices;
        {
            std::set<Cell> vertices_set;
            for (const Cell& edge : kspace.uLowerIncident(face))
                for (const Cell& vertex : kspace.uLowerIncident(edge))
                    if (calculus.containsCell(vertex))
                        vertices_set.insert(vertex);
            ASSERT( vertices_set.size() == 4 );
            std::copy(vertices_set.begin(), vertices_set.end(), std::back_inserter(vertices));
            std::swap(vertices[2], vertices[3]);
        }

        std::vector<RealPoint> vertex_positions;
        for (const Cell& vertex : vertices)
        {
            const Index vertex_index = calculus.getCellIndex(vertex);
            RealPoint position;
            for (int dim=0; dim<3; dim++)
                position[dim] = positions[vertex_index+dim*nvertices];
            vertex_positions.push_back(position);
        }

        for (int kk=0; kk<4; kk++)
        {
            const int kk_next = (kk+1)%4;
            align_energy += area*pow(face_normal.dot(vertex_positions[kk]-vertex_positions[kk_next]), 2);
        }
    }

    return std::make_tuple(position_energy, align_energy);
}

std::tuple<FlatVector, FlatVector, FlatVector, FlatVector>
approximateSurface(const Calculus& calculus, const FlatVector& normals, const ApproxParams& params)
{
    using DGtal::trace;
    using std::endl;
    using DGtal::PRIMAL;

    trace.beginBlock("approximating surface");

    typedef Calculus::Index Index;
    typedef std::vector<Triplet> Triplets;

    const KSpace& kspace = calculus.myKSpace;
    auto nfaces = calculus.kFormLength(2, PRIMAL);
    auto nedges = calculus.kFormLength(1, PRIMAL);
    auto nvertices = calculus.kFormLength(0, PRIMAL);
    trace.info() << "nfaces=" << nfaces << endl;
    trace.info() << "nedges=" << nedges << endl;
    trace.info() << "nvertices=" << nvertices << endl;

    trace.info() << "computing operator" << endl;

    FlatVector original_centers(3*nfaces);
    OperatorMatrix big_primal_hodge_2(3*nfaces, 3*nfaces);
    OperatorMatrix vertex_positions_to_face_centers(3*nfaces, 3*nvertices);
    OperatorMatrix align_normals(3*nvertices, 3*nvertices);
    OperatorMatrix fairness_operator(3*nvertices, 3*nvertices);
    {
        Triplets big_primal_hodge_2_triplets;
        Triplets vertex_positions_to_face_centers_triplets;
        Triplets align_normals_triplets;
        Triplets fairness_operator_triplets;
        for (const Calculus::ConstIterator::value_type& pair : calculus)
        {
            const Cell& face = pair.first;
            if (kspace.uDim(face) != 2) continue;
            const Index face_index = pair.second.index;
            const double face_area = pair.second.dual_size/pair.second.primal_size;
            ASSERT( face_area == 1 );

            { // primal hodge 2
                const double hodge_coeff = pair.second.dual_size/pair.second.primal_size;
                ASSERT( hodge_coeff == 1 );
                for (int dim=0; dim<3; dim++)
                    big_primal_hodge_2_triplets.push_back(Triplet(face_index+dim*nfaces, face_index+dim*nfaces, hodge_coeff));
            }

            std::set<Cell> face_vertices;
            for (const Cell& edge : kspace.uLowerIncident(face))
                for (const Cell& vertex : kspace.uLowerIncident(edge))
                    if (calculus.containsCell(vertex))
                        face_vertices.insert(vertex);
            ASSERT( face_vertices.size() == 4 );

#if !defined(NDEBUG)
            { // cell ordering check
                std::set<Cell>::const_iterator aa = face_vertices.begin();
                std::set<Cell>::const_iterator bb = face_vertices.begin();
                ASSERT( bb != face_vertices.end() );
                ++bb;
                ASSERT( bb != face_vertices.end() );
                std::vector<double> lengths;
                while (aa != face_vertices.end())
                {
                    const double length = (kspace.uKCoords(*aa)-kspace.uKCoords(*bb)).norm();
                    lengths.push_back(length);
                    ++aa;
                    ++bb;
                    if (bb == face_vertices.end()) bb= face_vertices.begin();
                }
                ASSERT( lengths == std::vector<double>({2,sqrt(8),2,sqrt(8)}) );
            }
#endif

            std::vector<Cell> ordered_face_vertices;
            std::copy(face_vertices.begin(), face_vertices.end(), std::back_inserter(ordered_face_vertices));
            std::swap(ordered_face_vertices[2], ordered_face_vertices[3]);

            { // align normal
                for (int normal_dim_ii=0; normal_dim_ii<3; normal_dim_ii++)
                    for (int normal_dim_jj=0; normal_dim_jj<3; normal_dim_jj++)
                    {
                        const double normal_coeff = face_area*normals(face_index+normal_dim_ii*nfaces)*normals(face_index+normal_dim_jj*nfaces);
                        //if (std::abs(normal_coeff) < 1e-3) continue;
                        for (int kk_vertex=0; kk_vertex<4; kk_vertex++)
                        {
                            const int kk_vertex_next = (kk_vertex+1)%4;
                            const int kk_vertex_prev = (kk_vertex+3)%4;
                            ASSERT( (kspace.uKCoords(ordered_face_vertices[kk_vertex])-kspace.uKCoords(ordered_face_vertices[kk_vertex_next])).norm() == 2 );
                            ASSERT( (kspace.uKCoords(ordered_face_vertices[kk_vertex])-kspace.uKCoords(ordered_face_vertices[kk_vertex_prev])).norm() == 2 );
                            const Index index_vertex = calculus.getCellIndex(ordered_face_vertices[kk_vertex]);
                            const Index index_vertex_next = calculus.getCellIndex(ordered_face_vertices[kk_vertex_next]);
                            const Index index_vertex_prev = calculus.getCellIndex(ordered_face_vertices[kk_vertex_prev]);
                            align_normals_triplets.push_back(Triplet(normal_dim_ii*nvertices+index_vertex, normal_dim_jj*nvertices+index_vertex, 2*normal_coeff));
                            align_normals_triplets.push_back(Triplet(normal_dim_ii*nvertices+index_vertex, normal_dim_jj*nvertices+index_vertex_prev, -normal_coeff));
                            align_normals_triplets.push_back(Triplet(normal_dim_ii*nvertices+index_vertex, normal_dim_jj*nvertices+index_vertex_next, -normal_coeff));
                        }
                    }
            }

            { // fairness operator
                for (int kk_vertex=0; kk_vertex<4; kk_vertex++)
                {
                    const Index kk_index = calculus.getCellIndex(ordered_face_vertices[kk_vertex]);
                    for (int ll_vertex=0; ll_vertex<4; ll_vertex++)
                    {
                        const Index ll_index = calculus.getCellIndex(ordered_face_vertices[ll_vertex]);
                        const double coeff = kk_vertex%2 == ll_vertex%2 ? face_area : -face_area;
                        for (int dim_kk=0; dim_kk<3; dim_kk++)
                            for (int dim_ll=0; dim_ll<3; dim_ll++)
                                fairness_operator_triplets.push_back(Triplet(kk_index+dim_kk*nvertices, ll_index+dim_kk*nvertices, coeff));
                    }
                }
            }

            { // original centers
                const CellEmbedder embedder(kspace);
                RealPoint original_center;
                for (const Cell& face_vertex : face_vertices)
                {
                    original_center += embedder(face_vertex);
                    const Index vertex_index = calculus.getCellIndex(face_vertex);
                    for (int dim=0; dim<3; dim++)
                        vertex_positions_to_face_centers_triplets.push_back(Triplet(face_index+dim*nfaces, vertex_index+dim*nvertices, 1./4));
                }
                original_center /= 4;

                for (int dim=0; dim<3; dim++) original_centers[face_index+dim*nfaces] = original_center[dim];
            }
        }
        big_primal_hodge_2.setFromTriplets(big_primal_hodge_2_triplets.begin(), big_primal_hodge_2_triplets.end());
        vertex_positions_to_face_centers.setFromTriplets(vertex_positions_to_face_centers_triplets.begin(), vertex_positions_to_face_centers_triplets.end());
        align_normals.setFromTriplets(align_normals_triplets.begin(), align_normals_triplets.end());
        fairness_operator.setFromTriplets(fairness_operator_triplets.begin(), fairness_operator_triplets.end());

    }

    FlatVector original_positions(3*nvertices);
    OperatorMatrix big_primal_hodge_0(3*nvertices, 3*nvertices);
    OperatorMatrix vertex_barycenter(3*nvertices, 3*nvertices);
    {
        const CellEmbedder embedder(kspace);
        Triplets big_primal_hodge_0_triplets;
        Triplets vertex_barycenter_triplets;
        for (const Calculus::ConstIterator::value_type& pair : calculus)
        {
            const Cell& vertex = pair.first;
            if (kspace.uDim(vertex) != 0) continue;
            const Index vertex_index = calculus.getCellIndex(vertex);

            { // primal hodge 0
                ASSERT( pair.second.primal_size == 1 );
                const double hodge_coeff = pair.second.dual_size/pair.second.primal_size;
                for (int dim=0; dim<3; dim++)
                    big_primal_hodge_0_triplets.push_back(Triplet(vertex_index+dim*nvertices, vertex_index+dim*nvertices, hodge_coeff));
            }

            { // original positions
                const RealPoint coords = embedder(vertex);
                for (int dim=0; dim<3; dim++) original_positions[pair.second.index+dim*nvertices] = coords[dim];
            }

            { // vertex barycenter
                std::set<Cell> neighbor_vertices;
                std::set<Cell> neighbor_edges;
                std::set<Cell> neighbor_faces;
                for (const Cell& edge : kspace.uUpperIncident(vertex))
                    if (calculus.containsCell(edge))
                    {
                        neighbor_edges.insert(edge);
                        for (const Cell& neighbor_vertex : kspace.uLowerIncident(edge))
                            if (neighbor_vertex != vertex && calculus.containsCell(neighbor_vertex))
                                neighbor_vertices.insert(neighbor_vertex);
                        for (const Cell& neighbor_face : kspace.uUpperIncident(edge))
                            if (calculus.containsCell(neighbor_face))
                                neighbor_faces.insert(neighbor_face);
                    }

                const bool on_border = (neighbor_edges.size()-1 == neighbor_faces.size()); // see dgtal_complex for enumeration of all cases

                if (!on_border)
                {
                    ASSERT( neighbor_vertices.size() >= 3 ); // >= 3 if surface has no border
                    ASSERT( neighbor_vertices.size() <= 6 );
                    if (neighbor_vertices.size() != pair.second.dual_size*4) trace.warning() << "non 2-manifold!!!!" << endl; // only detect manifold problem when surface has no border

                    const double coeff = 1./neighbor_vertices.size();
                    for (const Cell& neighbor_vertex : neighbor_vertices)
                    {
                        const Index neighbor_vertex_index = calculus.getCellIndex(neighbor_vertex);
                        for (int dim=0; dim<3; dim++)
                            vertex_barycenter_triplets.push_back(Triplet(vertex_index+dim*nvertices, neighbor_vertex_index+dim*nvertices, coeff));
                    }

                    for (int dim=0; dim<3; dim++)
                        vertex_barycenter_triplets.push_back(Triplet(vertex_index+dim*nvertices, vertex_index+dim*nvertices, -1));
                }
            }
        }
        big_primal_hodge_0.setFromTriplets(big_primal_hodge_0_triplets.begin(), big_primal_hodge_0_triplets.end());
        vertex_barycenter.setFromTriplets(vertex_barycenter_triplets.begin(), vertex_barycenter_triplets.end());
    }

    ASSERT( (vertex_positions_to_face_centers*original_positions-original_centers).array().maxCoeff() < 1e-5 );

    //typedef Backend::SolverSparseLU LinearSolver;
    typedef Backend::SolverSimplicialLDLT LinearSolver;
    const double lambda = params.regularization_position;
    const double alpha = params.regularization_center;
    const double gamma = params.align;
    const double beta = params.fairness;
    const double rau = params.barycenter;
    trace.info() << "lambda=" << lambda << " (vertex regularization)" << endl;
    trace.info() << "alpha=" << alpha << " (face center regularization)" << endl;
    trace.info() << "gamma=" << gamma << " (normal align)" << endl;
    trace.info() << "beta=" << beta << " (fairness)" << endl;
    trace.info() << "rau=" << rau << " (barycenter)" << endl;

    const OperatorMatrix ope =
    lambda*2*big_primal_hodge_0 +
    alpha*2*vertex_positions_to_face_centers.transpose()*big_primal_hodge_2*vertex_positions_to_face_centers +
    gamma*align_normals +
    beta*2*fairness_operator +
    rau*2*vertex_barycenter.transpose()*big_primal_hodge_0*vertex_barycenter;
    const FlatVector rht =
    lambda*2*big_primal_hodge_0*original_positions +
    alpha*2*vertex_positions_to_face_centers.transpose()*big_primal_hodge_2*original_centers;

    using namespace DGtal; // for solver status pretty print

    trace.info() << "factorizing " << (checkOperatorSymmetry(ope) ? "" : "non-") << "symmetric operator ";
    LinearSolver linear_solver;
    linear_solver.compute(ope);
    std::cout << linear_solver.info() << endl;

    if (linear_solver.info() == Eigen::Success)
    {
        trace.info() << "solving problem ";
        const FlatVector regularized_positions = linear_solver.solve(rht);
        std::cout << linear_solver.info() << endl;

        if (linear_solver.info() == Eigen::Success)
        {
            trace.endBlock();
            return std::make_tuple(original_positions, regularized_positions, original_centers, FlatVector(vertex_positions_to_face_centers*regularized_positions));
        }
    }

    trace.endBlock();
    return std::make_tuple(FlatVector(), FlatVector(), FlatVector(), FlatVector());
}

