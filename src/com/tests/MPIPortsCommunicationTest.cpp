#ifndef PRECICE_NO_MPI

#include "com/MPIPortsCommunication.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIPorts)

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcesses,
                     *testing::MinRanks(2)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>()
                     * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1}))
                     * boost::unit_test::label("MPI_Ports"))
{
  if (utils::Parallel::getCommunicatorSize() != 2)
    return;

  MPIPortsCommunication communication;

  std::string nameEven("even");
  std::string nameOdd("odd");

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    communication.acceptConnection(nameEven, nameOdd, 0, 1);
    int message = 1;
    communication.send(message, 0);
    communication.receive(message, 0);
    BOOST_TEST(message == 2);
    communication.closeConnection();
    break;
  }
  case 1: {
    communication.requestConnection(nameEven, nameOdd, 0, 1);
    int message = -1;
    communication.receive(message, 0);
    BOOST_TEST(message == 1);
    message = 2;
    communication.send(message, 0);
    communication.closeConnection();
    break;
  }
  }
}

BOOST_AUTO_TEST_SUITE_END() // MPIPortsCommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI
