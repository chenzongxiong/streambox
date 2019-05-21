#ifndef APPLIEDPTRANSFORM_H_
#define APPLIEDPTRANSFORM_H_
#include <assert.h>
#include <string>
using namespace std;
class PValue;
//template <class Derived>
class PTransform;
////////////////////////////////////////////////////////////////////////

class AppliedPTransform {

private:
	string _fullName;

	// note we can't simply make copies as they could point to subclasses
	PValue* _input;
	PValue* _output;

	PTransform *_transform;

	AppliedPTransform(string fullName, PValue* input, PValue* output, PTransform* transform) :
		_fullName(fullName), _input(input), _output(output), _transform(transform) { }

public:
	// pass in references. return a reference
	static AppliedPTransform* of(string fullName, PValue* input, PValue* output, PTransform* transform) {
		return new AppliedPTransform(fullName, input, output, transform);
	}

	string getFullName() { return _fullName; }

	PValue* getInput() { return _input; }

	PValue* getOutput() { return _output; }

	PTransform* getTransform() { return _transform; }

	int hashCode() { assert(0); return 0; } // TODO

	bool operator==(AppliedPTransform& other) {
		// TODO: better check
		return ((this->getFullName() == other.getFullName())
				&& (this->getInput() == other.getInput())
				&& (this->getOutput() == other.getOutput())
				&& (this->getTransform() == other.getTransform()));
	}
};

#endif /* APPLIEDPTRANSFORM_H_ */
