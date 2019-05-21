#ifndef WINDOWINGSTRATEGY_H_
#define WINDOWINGSTRATEGY_H_
#include "BoundedWindow.h"
class WindowingStrategy : public BoundedWindow {
public:
	enum AccumulationMode {
	    DISCARDING_FIRED_PANES,
	    ACCUMULATING_FIRED_PANES
	};

private:
	WindowingStrategy() {
		// TODO
	}
};

#endif /* WINDOWINGSTRATEGY_H_ */
