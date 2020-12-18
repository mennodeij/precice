#pragma once

#include <map>
#include <set>
#include <stddef.h>
#include <string>
#include <vector>
#include "action/Action.hpp"
#include "boost/noncopyable.hpp"
#include "com/Communication.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "io/Constants.hpp"
#include "logging/Logger.hpp"
#include "m2n/BoundM2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/DataContext.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "utils/MultiLock.hpp"

namespace precice {
namespace config {
class SolverInterfaceConfiguration;
}
} // namespace precice

// Forward declaration to friend the boost test struct
namespace PreciceTests {
namespace Serial {
struct TestConfigurationPeano;
struct TestConfigurationComsol;
} // namespace Serial
} // namespace PreciceTests

namespace precice {
namespace cplscheme {
class CouplingSchemeConfiguration;
} // namespace cplscheme
namespace mesh {
class Mesh;
} // namespace mesh

namespace impl {

/// Implementation of solver interface.
class SolverInterfaceImpl {
public:
  /**
   * @brief Constructor.
   *
   * A solver that wants to use the SolverInterfaceImpl must instatiate an object
   * of this class. The object has to be configured by one of the configure
   * methods before it has a reasonable state and can be used.
   *
   * @param[in] configurationFileName Name (with path) of the xml config file.
   * @param[in] participantName Name of the participant using the interface. Has to
   *                            match the name given for a participant in the
   *                            xml configuration file.
   */
  SolverInterfaceImpl(
      std::string        participantName,
      const std::string &configurationFileName,
      int                accessorProcessRank,
      int                accessorCommunicatorSize);

  /// Deleted copy constructor
  SolverInterfaceImpl(SolverInterfaceImpl const &) = delete;

  /// Deleted copy assignment
  SolverInterfaceImpl &operator=(SolverInterfaceImpl const &) = delete;

  /// Deleted move constructor
  SolverInterfaceImpl(SolverInterfaceImpl &&) = delete;

  /// Deleted move assignment
  SolverInterfaceImpl &operator=(SolverInterfaceImpl &&) = delete;

  /**
   * @brief Constructor with support for custom MPI_COMM_WORLD.
   *
   * A solver that wants to use the SolverInterfaceImpl must instatiate an object
   * of this class. The object has to be configured by one of the configure
   * methods before it has a reasonable state and can be used.
   *
   * Use the parameter communicator to specify a custom global MPI communicator.
   * Pass a null pointer to signal preCICE to use MPI_COMM_WORLD.
   *
   * @param[in] participantName Name of the participant using the interface. Has to
   *                            match the name given for a participant in the
   *                            xml configuration file.
   * @param[in] configurationFileName Name (with path) of the xml config file.
   * @param[in] communicator    A pointer to the MPI_Comm to use.
   */
  SolverInterfaceImpl(
      std::string        participantName,
      const std::string &configurationFileName,
      int                accessorProcessRank,
      int                accessorCommunicatorSize,
      void *             communicator);

  /** Ensures that finalize() has been called.
   *
   * @see finalize()
   */
  ~SolverInterfaceImpl();

  /**
   * @brief Initializes all coupling data and starts (coupled) simulation.
   *
   * - Initiates MPI communication, if not done yet.
   * - Sets up a connection to the other participants of the coupled simulation.
   * - Creates all meshes, solver meshes need to be submitted before.
   * - Receives first coupling data, when the solver is not starting the
   *   coupled simulation.
   * - Determines length of the first timestep to be computed.
   *
   * @return Maximum length of first timestep to be computed by the solver.
   */
  double initialize();

  /**
   * @brief Sets the sofar written data as initial value for the coupling scheme.
   *
   * Erases the written data afterwards.
   */
  void initializeData();

  /**
   * @brief Exchanges coupling data and advances coupling state.
   *
   * - Sends and resets coupling data written by solver to coupling partners.
   * - Receives coupling data read by solver.
   * - Computes and applied data mappings.
   * - Computes acceleration of coupling data.
   * - Exchanges and computes information regarding the state of the coupled
   *   simulation.
   *
   * @param[in] computedTimestepLength Length of timestep computed by solver.
   * @return Maximum length of next timestep to be computed by solver.
   */
  double advance(double computedTimestepLength);

  /**
   * @brief Finalizes the coupled simulation.
   *
   * If initialize() has been called:
   *
   * - Synchronizes with remote partiticipants
   * - handles final exports
   * - cleans up general state
   *
   * Always:
   *
   * - flushes and finalizes Events
   * - finalizes managed PETSc
   * - finalizes managed MPI
   *
   * @post Closes MasterSlave communication
   * @post Finalized managed PETSc
   * @post Finalized managed MPI
   *
   * @warning
   * Finalize is not the inverse of initialize().
   * Finalize has to be called after construction.
   */
  void finalize();

  /**
   * @brief Returns the number of spatial dimensions for the coupling.
   *
   * The number of dimensions is fixed on configuration.
   */
  int getDimensions() const;

  /**
   * @brief Returns true, if the coupled simulation is still ongoing.
   *
   * The information to decide about the continuation of the coupled simulation
   * is retreived in the function initializeCoupling and updated in the
   * function exchangeData.
   */
  bool isCouplingOngoing() const;

  /**
   * @brief Returns true, if new data to be read is available.
   */
  bool isReadDataAvailable() const;

  /**
   * @brief Returns true, if new data has to be written.
   */
  bool isWriteDataRequired(double computedTimestepLength) const;

  /**
   * @brief Returns true, if a time window is completed.
   */
  bool isTimeWindowComplete() const;

  /**
   * @brief Returns whether the solver has to evaluate the surrogate model representation
   *        It does not automatically imply, that the solver does not have to evaluate the
   *        fine model representation
   */
  bool hasToEvaluateSurrogateModel() const;

  /**
   * @brief Returns whether the solver has to evaluate the fine model representation
   *        It does not automatically imply, that the solver does not have to evaluate the
   *        surrogate model representation
   */
  bool hasToEvaluateFineModel() const;

  /**
   * @brief Returns true, if provided name of action is required.
   *
   * Some features of preCICE require a solver to perform specific actions, in
   * order to be in valid state for a coupled simulation. A solver is made
   * eligible to use those features, by querying for the required actions,
   * performing them on demand, and calling markActionFulfilled() to signalize
   * preCICE the correct behavior of the solver.
   */
  bool isActionRequired(const std::string &action) const;

  /**
   * @brief Tells preCICE that a required action has been fulfilled by a solver.
   *
   * For more details see method requireAction().
   */
  void markActionFulfilled(const std::string &action);

  /// Returns true, if the mesh with given name is used.
  bool hasMesh(const std::string &meshName) const;

  /**
   * @brief Returns the ID belonging to the mesh with given name.
   *
   * The existing names are determined from the configuration.
   */
  int getMeshID(const std::string &meshName) const;

  /// Returns all mesh IDs (besides sub-ids).
  std::set<int> getMeshIDs() const;

  /// Returns true, if the data with given name is used in the given mesh.
  bool hasData(const std::string &dataName, int meshID) const;

  /// Returns data id corresponding to the given name (from configuration) and mesh.
  int getDataID(const std::string &dataName, int meshID) const;

  /**
   * @brief Checks if the gradient with given data name is used by a solver and mesh.
   *
   * @param[in] dataName the name of the data
   * @param[in] meshID the id of the associated mesh
   * @returns whether the gradient is used.
   */
  bool hasGradient(const std::string &dataName, int meshID) const;

  /// Returns the number of nodes of a mesh.
  int getMeshVertexSize(int meshID) const;

  /**
   * @brief Resets mesh with given ID.
   *
   * Has to be called, everytime the positions for data to be mapped
   * changes. Only has an effect, if the mapping used is non-stationary and
   * non-incremental.
   */
  void resetMesh(int meshID);

  /**
   * @brief Set the position of a solver mesh vertex.
   *
   * @return Vertex ID to be used when setting an edge.
   */
  int setMeshVertex(
      int           meshID,
      const double *position);

  /**
   * @brief Sets several spatial positions for a mesh.
   *
   * @param[out] ids IDs for data from given positions.
   */
  void setMeshVertices(
      int           meshID,
      int           size,
      const double *positions,
      int *         ids);

  /**
   * @brief Gets spatial positions of vertices for given IDs.
   *
   * @param[in] ids IDs obtained when setting write positions.
   * @param[out] positions Positions corresponding to IDs.
   */
  void getMeshVertices(
      int        meshID,
      size_t     size,
      const int *ids,
      double *   positions) const;

  /**
   * @brief Gets vertex data ids from positions.
   *
   * @param[in] size Number of positions, ids.
   * @param[in] positions Positions (x,y,z,x,y,z,...) to find ids for.
   * @param[out] ids IDs corresponding to positions.
   */
  void getMeshVertexIDsFromPositions(
      int           meshID,
      size_t        size,
      const double *positions,
      int *         ids) const;

  /**
   * @brief Set an edge of a solver mesh.
   *
   * @return Index of the edge to be used when setting a triangle.
   */
  int setMeshEdge(
      int meshID,
      int firstVertexID,
      int secondVertexID);

  /// Set a triangle of a solver mesh.
  void setMeshTriangle(
      int meshID,
      int firstEdgeID,
      int secondEdgeID,
      int thirdEdgeID);

  /// Sets a triangle and creates/sets edges automatically of a solver mesh.
  void setMeshTriangleWithEdges(
      int meshID,
      int firstVertexID,
      int secondVertexID,
      int thirdVertexID);

  /// Set a quadrangle of a solver mesh.
  void setMeshQuad(
      int meshID,
      int firstEdgeID,
      int secondEdgeID,
      int thirdEdgeID,
      int fourthEdgeID);

  /// Sets a quadrangle and creates/sets edges automatically of a solver mesh.
  void setMeshQuadWithEdges(
      int meshID,
      int firstVertexID,
      int secondVertexID,
      int thirdVertexID,
      int fourthVertexID);

  /**
   * @brief Computes and maps all write data mapped from mesh with given ID.
   *
   * Is automatically called in advance, if not called manually before.
   */
  void mapWriteDataFrom(int fromMeshID);

  /// Computes and maps all read data mapped to mesh with given ID.
  void mapReadDataTo(int toMeshID);

  /**
   * @brief Writes vector data values given as block.
   *
   * The block must contain the vector values in the following form:
   * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
   * the number of vector values. In 2D, the z-components are removed.
   *
   * @param[in] fromDataID ID of the data to be written.
   * @param[in] size Number of valueIndices, and number of values * dimensions.
   * @param[in] values Values of the data to be written.
   */
  void writeBlockVectorData(
      int           fromDataID,
      int           size,
      const int *   valueIndices,
      const double *values);

  /**
   * @brief Write vectorial data to the interface mesh
   *
   * The exact mapping and communication must be specified in XYZ.
   *
   * @param[in] fromDataID ID of the data to be written, e.g. 1 = forces
   * @param[in] dataPosition Position (coordinate, e.g.) of data to be written
   * @param[in] dataValue Value of the data to be written
   */
  void writeVectorData(
      int           fromDataID,
      int           valueIndex,
      const double *value);

  /**
   * @brief Writes gradient given as block.
   *
   * This function writes values of specified vertices to a dataID.
   * Values are provided as a block of continuous memory.
   * valueIndices contains the indices of the vertices
   *
   * For each vertex, the gradient is an m-by-n matrix, where the number
   * of rows represent the value dimension, and the number of columns
   * the spatial dimension. For example, the gradient of a scalar in 3D is
   * given in a 1x3 matrix, which when multiplied by the location difference
   * column vector (i.e. a 3x1 matrix) gives the gradient contribution as a
   * scalar. Usually, for vector fields the value dimension and spatial dimension are
   * equal, but this is not required.
   *
   * The format in 2 and 3 spatial dimensions is as follows:
   *   Scalar: Field on each vertex is (s)
   *      ds0/dx, ds0/dy, ds1/dx, ds1/dy, ... , dsN/dx, dsN/dy (2D)
   *      ds0/dx, ds0/dy, ds0/dz, ds1/dx, ds1/dy, ds1/dz, ... , dsN/dx, dsN/dy, dsN/dz (3D)
   *
   *  Vector2: Field on each vertex is (u,v)
   *      du0/dx, du0/dy, du1/dx, du1/dy, ..., duN/dx, duN/dy,
   *      dv0/dx, dv0/dy, dv1/dx, dv1/dy, ..., dvN/dx, dvN/dy (2D)
   *
   *      du0/dx, du0/dy, du0/dz, du1/dx, du1/dy, du1/dz, ..., duN/dx, duN/dy, duN/dz,
   *      dv0/dx, dv0/dy, dv0/dz, dv1/dx, dv1/dy, dv1/dz, ..., dvN/dx, dvN/dy, dvN/dz (3D)
   *
   *  Vector3: Field on each vertex is (u,v,w)
   *      du0/dx, du0/dy, du1/dx, du1/dy, ..., duN/dx, duN/dy,
   *      dv0/dx, dv0/dy, dv1/dx, dv1/dy, ..., dvN/dx, dvN/dy,
   *      dw0/dx, dw0/dy, dw1/dx, dw1/dy, ..., dwN/dx, dwN/dy (2D)
   *
   *      du0/dx, du0/dy, du0/dz, du1/dx, du1/dy, du1/dz, ..., duN/dx, duN/dy, duN/dz,
   *      dv0/dx, dv0/dy, dv0/dz, dv1/dx, dv1/dy, dv1/dz, ..., dvN/dx, dvN/dy, dvN/dz,
   *      dw0/dx, dw0/dy, dw0/dz, dw1/dx, dw1/dy, dw1/dz, ..., dwN/dx, dwN/dy, dwN/dz (3D)
   *
   * @param[in] forDataID data ID to write gradients for.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[in] values pointer to the vector values.
   *
   * @pre count of available elements at values matches the configured dimension * size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeBlockGradient(
      int           forDataID,
      int           size,
      const int*    valueIndices,
      const double* gradientValues);

  /**
   * @brief Writes gradient to a vertex
   *
   * The format in 2 and 3 spatial dimensions is as follows:
   *   Scalar: Field on the vertex is (s)
   *      ds/dx, ds/dy        (2D)
   *      ds/dx, ds/dy, ds/dz (3D)
   *
   *  Vector2: Field on each vertex is (u,v)
   *      du/dx, du/dy,
   *      dv/dx, dv/dy, (2D)
   *
   *      du/dx, du/dy, du/dz,
   *      dv/dx, dv/dy, dv/dz(3D)
   *
   *  Vector3: Field on each vertex is (u,v,w)
   *      du/dx, du/dy,
   *      dv/dx, dv/dy,
   *      dw/dx, dw/dy (2D)
   *
   *      du/dx, du/dy, du/dz,
   *      dv/dx, dv/dy, dv/dz,
   *      dw/dx, dw/dy, dw/dz (3D)
   *
   * @param[in] forDataID data ID to write gradients for.
   * @param[in] valueIndex Index of the vertex.
   * @param[in] value pointer to the gradient value.
   */
  void writeGradient(
      int           forDataID,
      int           valueIndex,
      const double *value);

  /**
   * @brief Writes scalar data values given as block.
   *
   * @param fromDataID [IN] ID of the data to be written.
   * @param size [IN] Number of valueIndices, and number of values.
   * @param values [IN] Values of the data to be written.
   */
  void writeBlockScalarData(
      int           fromDataID,
      int           size,
      const int *   valueIndices,
      const double *values);

  /**
   * @brief Write scalar data to the interface mesh
   *
   * The exact mapping and communication must be specified in XYZ.
   *
   * @param fromDataID       [IN] ID of the data to be written (2 = temperature, e.g.)
   * @param dataPosition [IN] Position (coordinate, e.g.) of data to be written
   * @param dataValue    [IN] Value of the data to be written
   */
  void writeScalarData(
      int    fromDataID,
      int    valueIndex,
      double value);

  /**
   * @brief Reads vector data values given as block.
   *
   * The block contains the vector values in the following form:
   * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
   * the number of vector values. In 2D, the z-components are removed.
   *
   * @param toDataID [IN] ID of the data to be read.
   * @param size [IN] Number of indices, and number of values * dimensions.
   * @param valueIndices [IN] Indices (from setReadPosition()) of data values.
   * @param values [IN] Values of the data to be read.
   */
  void readBlockVectorData(
      int        toDataID,
      int        size,
      const int *valueIndices,
      double *   values) const;

  /**
   * @brief Reads vector data from the coupling mesh.
   *
   * @param[in] toDataID ID of the data to be read, e.g. 1 = forces
   * @param[in] dataPosition Position (coordinate, e.g.) of data to be read
   * @param[out] dataValue Read data value
   */
  void readVectorData(
      int     toDataID,
      int     valueIndex,
      double *value) const;

  /**
   * @brief Reads scalar data values given as block.
   *
   * @param[in] toDataID ID of the data to be written.
   * @param[in] size Number of valueIndices, and number of values.
   * @param[in] values Values of the data to be written.
   */
  void readBlockScalarData(
      int        toDataID,
      int        size,
      const int *valueIndices,
      double *   values) const;

  /**
   * @brief Read scalar data from the interface mesh.
   *
   * The exact mapping and communication must be specified in XYZ.
   *
   * @param[in] toDataID     ID of the data to be read, e.g. 2 = temperatures
   * @param[in] dataPosition Position (coordinate, e.g.) of data to be read
   * @param[in] dataValue    Read data value
   */
  void readScalarData(
      int     toDataID,
      int     valueIndex,
      double &value) const;

  /**
   * @brief Sets the location for all output of preCICE.
   *
   * If done after configuration, this overwrites the output location specified
   * in the configuration.
   */
  //  void setExportLocation (
  //    const std::string& location,
  //    int                exportType = io::constants::exportAll() );

  /**
   * @brief Writes a mesh to vtk file.
   *
   * The plotting path has to be specified in the configuration of the
   * accessing participant.
   *
   * @param[in] filenameSuffix Suffix of all plotted files
   */
  void exportMesh(
      const std::string &filenameSuffix,
      int                exportType = io::constants::exportAll()) const;

  /**
   * @brief Scales data values according to configuration.
   *
   * Currently, the only scaling supported is a division of the data values
   * through the surface area belonging to its "support". This allows to come
   * from forces to stresses, e.g..
   */
  //  void scaleReadData ()

  /// Returns the name of the accessor
  std::string getAccessorName() const
  {
    return _accessorName;
  }

  /// Returns the rank of the accessor
  int getAccessorProcessRank() const
  {
    return _accessorProcessRank;
  }

  /// Returns the size of the accessors communicator
  int getAccessorCommunicatorSize() const
  {
    return _accessorCommunicatorSize;
  }

  /// Allows to access a registered mesh
  const mesh::Mesh &mesh(const std::string &meshName) const;

private:
  mutable logging::Logger _log{"impl::SolverInterfaceImpl"};

  std::string _accessorName;

  int _accessorProcessRank;

  int _accessorCommunicatorSize;

  impl::PtrParticipant _accessor;

  /// Spatial dimensions of problem.
  int _dimensions = 0;

  utils::MultiLock<int> _meshLock;

  /// mesh name to mesh ID mapping.
  std::map<std::string, int> _meshIDs;

  /// dataIDs referenced by meshID and data name
  std::map<int, std::map<std::string, int>> _dataIDs;

  /// gradientIDs referenced by meshID and data name
  std::map<int, std::map<std::string, int>> _gradientIDs;

  std::map<std::string, m2n::BoundM2N> _m2ns;

  /// Holds information about solvers participating in the coupled simulation.
  std::vector<impl::PtrParticipant> _participants;

  cplscheme::PtrCouplingScheme _couplingScheme;

  /// Represents the various states a SolverInterface can be in.
  enum struct State {
    Constructed, // Initial state of SolverInterface
    Initialized, // SolverInterface.initialize() triggers transition from State::Constructed to State::Initialized; mandatory
    Finalized    // SolverInterface.finalize() triggers transition form State::Initialized or State::InitializedData to State::Finalized; mandatory
  };

  // SolverInterface.initializeData() triggers transition from false to true.
  bool _hasInitializedData = false;

  /// The current State of the solverinterface
  State _state{State::Constructed};

  /// Counts calls to advance for plotting.
  long int _numberAdvanceCalls = 0;

  /**
   * @brief Configures the coupling interface from the given xml file.
   *
   * Only after the configuration a reasonable state of a SolverInterfaceImpl
   * object is achieved.
   *
   * @param configurationFileName [IN] Name (with path) of the xml config. file.
   */
  void configure(const std::string &configurationFileName);

  /**
   * @brief Configures the coupling interface with a prepared configuration.
   *
   * Can be used to configure the SolverInterfaceImpl without xml file. Requires
   * to manually setup the configuration object.
   */
  void configure(const config::SolverInterfaceConfiguration &configuration);

  void configureM2Ns(const m2n::M2NConfiguration::SharedPointer &config);

  /// Exports meshes with data and watch point data.
  void handleExports();

  /**
   * @brief Adds exchanged data ids related to accessor to the coupling scheme.
   *
   * From the configuration it is only known which data a participant writes,
   * which data he receives, and which data is exchanged within the coupling
   * scheme. The data actually being sent and received in the coupling scheme,
   * is defined by the intersection of data being written or read by the
   * accessor and the data exchanged in the coupling scheme, respectively.
   *
   * Prerequesits:
   * - _accessor is valid
   * - _couplingScheme is valid
   */
  void addDataToCouplingScheme(
      const cplscheme::CouplingSchemeConfiguration &config);

  /**
   * @brief Configures _couplingScheme to be ready to use.
   *
   * @pre _couplingScheme holds a pointer to a valid coupling scheme object
   * @pre all meshes are created
   * @pre all data is added to _data
   * @pre _accessor points to the accessor
   */
  void addMeshesToCouplingScheme();

  /// Returns true, if the accessor uses the mesh with given name.
  bool isUsingMesh(const std::string &meshName) const;

  /// Determines participants providing meshes to other participants.
  void configurePartitions(
      const m2n::M2NConfiguration::SharedPointer &m2nConfig);

  /// Communicate bounding boxes and look for overlaps
  void compareBoundingBoxes();

  /// Communicate meshes and create partitions
  void computePartitions();

  /// Computes, performs, and resets all suitable write mappings.
  void mapWrittenData();

  /// Computes, performs, and resets all suitable read mappings.
  void mapReadData();

  /**
   * @brief Performs all data actions with given timing.
   *
   * @param[in] dt Last timestep length computed by solver.
   * @param[in] partFullDt Part of current full timestep computed by solver.
   * @param[out] fullDt Length of current full timestep.
   */
  void performDataActions(
      const std::set<action::Action::Timing> &timings,
      double                                  time,
      double                                  dt,
      double                                  partFullDt,
      double                                  fullDt);

  /// Resets written data, displacements and mesh neighbors to export.
  void resetWrittenData();

  /// Determines participant accessing this interface from the configuration.
  impl::PtrParticipant determineAccessingParticipant(
      const config::SolverInterfaceConfiguration &config);

  int markedSkip() const
  {
    return 0;
  }
  int markedQueryDirectly() const
  {
    return 1;
  }

  /// Initializes communication between master and slaves.
  void initializeMasterSlaveCommunication();

  /// Syncs the timestep between slaves and master (all timesteps should be the same!)
  void syncTimestep(double computedTimestepLength);

  /// To allow white box tests.
  friend struct PreciceTests::Serial::TestConfigurationPeano;
  friend struct PreciceTests::Serial::TestConfigurationComsol;
};

} // namespace impl
} // namespace precice
