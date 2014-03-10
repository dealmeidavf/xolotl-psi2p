#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DummyEventCounter.h>
#include <string>

using namespace std;
using namespace xolotlPerf;

/**
 * This suite is responsible for testing the DummyEventCounter.
 */
BOOST_AUTO_TEST_SUITE (DummyEventCounter_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

	DummyEventCounter tester("test");

	BOOST_REQUIRE_EQUAL("test", tester.getName());
}

BOOST_AUTO_TEST_CASE(checkInitialValue) {

	DummyEventCounter tester("test");

	BOOST_REQUIRE_EQUAL(0, tester.getValue());

}

BOOST_AUTO_TEST_CASE(checkCounting) {

	DummyEventCounter tester("test");

	long count = 3;

	for(int i = 0; i < 3; i++){
		tester.increment();
	}

	BOOST_REQUIRE_EQUAL(0, tester.getValue());

}


BOOST_AUTO_TEST_SUITE_END()

