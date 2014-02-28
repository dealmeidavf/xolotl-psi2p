#ifndef DUMMYTIMER_H
#define DUMMYTIMER_H

#include "ITimer.h"
#include "Identifiable.h"

using namespace std;

namespace xolotlPerf{

/**
 * The DummyTimer class is instantiated by the DummerHandlerRegistry class
 * and realizes the DummyTimer interface.
 */
class DummyTimer : public ITimer, public xolotlCore::Identifiable
{
private:

	/**
	 * The default constructor is declared as private since Timers
	 *  must be initialized with a name.
	 */
    DummyTimer()
      : xolotlCore::Identifiable("unused")
    { }

public:

	/**
	 * DummyTimer constructor that takes the argument timerName
	 * to distinguish specific DummyTimer.
	 *
	 * @param timerName The DummyTimer's name
	 */
	DummyTimer(std::string timerName)
      : xolotlCore::Identifiable("unused")
    { }

	/**
	 * The destructor.
	 */
	virtual ~DummyTimer() { }

    /**
     * This operations starts the ITimer.
     */
	virtual void start();

    /**
     * This operation stops the ITimer.
     */
	virtual void stop();

    /**
     * This operation returns the value of the DummyTimer.
     */
    virtual double getValue() const;

	/**
	 * This operation returns the units of the GPTLTimer.
	 */
    virtual long getUnits() const;

};  //end class DummyTimer

}  //end namespace xolotlPerf

#endif
