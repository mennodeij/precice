#include "NearestNeighborGradientMapping.hpp"

#include <Eigen/Core>
#include <boost/container/flat_set.hpp>
#include <functional>
#include <memory>

#include <boost/version.hpp>
#if BOOST_VERSION < 106600
#include <boost/function_output_iterator.hpp>
#else
#include <boost/iterator/function_output_iterator.hpp>
#endif

#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Gradient.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/RTree.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Event.hpp"
#include "utils/Statistics.hpp"
#include "utils/assertion.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

NearestNeighborGradientMapping::NearestNeighborGradientMapping(
    Constraint constraint,
    int        dimensions)
    : Mapping(constraint, dimensions)
{
  setInputRequirement(Mapping::MeshRequirement::VERTEX);
  setOutputRequirement(Mapping::MeshRequirement::VERTEX);
}

void NearestNeighborGradientMapping::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  const std::string     baseEvent = "map.nn.computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::utils::Event e(baseEvent, precice::syncMode);

  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
    precice::utils::Event e2(baseEvent + ".getIndexOnVertices", precice::syncMode);
    auto                  rtree = mesh::rtree::getVertexRTree(input());
    e2.stop();
    size_t verticesSize = output()->vertices().size();
    _vertexIndices.resize(verticesSize);
    utils::statistics::DistanceAccumulator distanceStatistics;
    const mesh::Mesh::VertexContainer &    outputVertices = output()->vertices();
    for (size_t i = 0; i < verticesSize; i++) {
      const Eigen::VectorXd &coords = outputVertices[i].getCoords();
      // Search for the output vertex inside the input mesh and add index to _vertexIndices
      rtree->query(boost::geometry::index::nearest(coords, 1),
                   boost::make_function_output_iterator([&](size_t const &val) {
                     const auto &match = input()->vertices()[val];
                     _vertexIndices[i] = match.getID();
                     distanceStatistics(bg::distance(match, coords));
                   }));
    }
    if (distanceStatistics.empty()) {
      PRECICE_INFO("Mapping distance not available due to empty partition.");
    } else {
      PRECICE_INFO("Mapping distance " << distanceStatistics);
    }
  } else {
    PRECICE_ASSERT(false, "Conservative mapping not available for nearest-neighbor-gradient mapping");
    //PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    //PRECICE_DEBUG("Compute conservative mapping");
    //precice::utils::Event e2(baseEvent + ".getIndexOnVertices", precice::syncMode);
    //auto                  rtree = mesh::rtree::getVertexRTree(output());
    //e2.stop();
    //size_t verticesSize = input()->vertices().size();
    //_vertexIndices.resize(verticesSize);
    //utils::statistics::DistanceAccumulator distanceStatistics;
    //const mesh::Mesh::VertexContainer &    inputVertices = input()->vertices();
    //for (size_t i = 0; i < verticesSize; i++) {
    //  const Eigen::VectorXd &coords = inputVertices[i].getCoords();
    //  // Search for the input vertex inside the output mesh and add index to _vertexIndices
    //  rtree->query(boost::geometry::index::nearest(coords, 1),
    //               boost::make_function_output_iterator([&](size_t const &val) {
    //                 const auto &match = output()->vertices()[val];
    //                 _vertexIndices[i] = match.getID();
    //                 distanceStatistics(bg::distance(match, coords));
    //               }));
    //}
    //if (distanceStatistics.empty()) {
    //  PRECICE_INFO("Mapping distance not available due to empty partition.");
    //} else {
    //  PRECICE_INFO("Mapping distance " << distanceStatistics);
    //}
  }
  _hasComputedMapping = true;
}

bool NearestNeighborGradientMapping::hasComputedMapping() const
{
  PRECICE_TRACE(_hasComputedMapping);
  return _hasComputedMapping;
}

void NearestNeighborGradientMapping::clear()
{
  PRECICE_TRACE();
  _vertexIndices.clear();
  _hasComputedMapping = false;
  if (getConstraint() == CONSISTENT) {
    mesh::rtree::clear(*input());
  } else {
    mesh::rtree::clear(*output());
  }
}

void NearestNeighborGradientMapping::map(
    int inputDataID,
    int outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);

  auto &data = input()->data(inputDataID);
  auto &gradient = input()->gradient(data);
  // note: no need for gradient != nullptr assert, the Mesh::gradient(const PtrData) does a sanity check

  precice::utils::Event e("map.nn.mapData.From" + input()->getName() + "To" + output()->getName() + " using gradient", precice::syncMode);

  const Eigen::VectorXd &inputValues   = data->values();
  const Eigen::MatrixXd &inputGradient = gradient->values();
  Eigen::VectorXd       &outputValues  = output()->data(outputDataID)->values();

  int valueDimensions = input()->data(inputDataID)->getDimensions();
  int spaceDimensions = input()->getDimensions();
  PRECICE_ASSERT(valueDimensions == output()->data(outputDataID)->getDimensions(),
                 valueDimensions, output()->data(outputDataID)->getDimensions());
  PRECICE_ASSERT(inputValues.size() / valueDimensions == (int) input()->vertices().size(),
                 inputValues.size(), valueDimensions, input()->vertices().size());
  PRECICE_ASSERT(outputValues.size() / valueDimensions == (int) output()->vertices().size(),
                 outputValues.size(), valueDimensions, output()->vertices().size());
  PRECICE_ASSERT(inputGradient.cols() / spaceDimensions == (int) input()->vertices().size(),
                 inputGradient.cols(), spaceDimensions, input()->vertices().size());
  PRECICE_ASSERT(inputGradient.rows() == valueDimensions,
                 inputGradient.rows(), valueDimensions);


  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Map consistent");
    size_t const outSize = output()->vertices().size();
    for (size_t i = 0; i < outSize; i++) {
      int inputIndex = _vertexIndices[i];// * valueDimensions;

      // get a view on the required gradient tensor
      const Eigen::MatrixXd &grad  = inputGradient.block(0,i*spaceDimensions,valueDimensions,spaceDimensions);
      PRECICE_DEBUG("Gradient matrix at inputIndex " << inputIndex << " = " << grad);
      // compute the distance between the required point and the found point
      const Eigen::VectorXd &delta = output()->vertices()[inputIndex].getCoords() - input()->vertices()[i].getCoords();

      PRECICE_ASSERT(delta.size() == spaceDimensions, delta.size(), spaceDimensions);

      // get the input value 
      const Eigen::VectorXd &value = inputValues.block(inputIndex*valueDimensions, 0, valueDimensions, 1);

      PRECICE_ASSERT(value.size() == valueDimensions, value.size(), valueDimensions);

      // compute gradient x vector product
      Eigen::VectorXd product = grad*delta;

      PRECICE_ASSERT(product.size() == valueDimensions, value.size(), valueDimensions);

      // compute the output value
      Eigen::VectorXd out = value + product;

      // store output value
      outputValues.block(i*valueDimensions, 0, valueDimensions, 1) = out;
    }
  } else {
    // does a CONSERVATIVE mapping even make sense, there is nothing to conserve?
    PRECICE_ASSERT(false, "CONSERVATIVE gradient mapping is not implemented");
    //PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    //PRECICE_DEBUG("Map conservative");
    //size_t const inSize = input()->vertices().size();
    //for (size_t i = 0; i < inSize; i++) {
    //  int const outputIndex = _vertexIndices[i] * valueDimensions;
    //  for (int dim = 0; dim < valueDimensions; dim++) {
    //    outputValues(outputIndex + dim) += inputValues((i * valueDimensions) + dim);
    //  }
    //}
  }
}

void NearestNeighborGradientMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();
  precice::utils::Event e("map.nn.tagMeshFirstRound.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  computeMapping();

  // Lookup table of all indices used in the mapping
  const boost::container::flat_set<int> indexSet(_vertexIndices.begin(), _vertexIndices.end());

  if (getConstraint() == CONSISTENT) {
    for (mesh::Vertex &v : input()->vertices()) {
      if (indexSet.count(v.getID()) != 0)
        v.tag();
    }
  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    for (mesh::Vertex &v : output()->vertices()) {
      if (indexSet.count(v.getID()) != 0)
        v.tag();
    }
  }

  clear();
}

void NearestNeighborGradientMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // for NN mapping no operation needed here
}

} // namespace mapping
} // namespace precice
