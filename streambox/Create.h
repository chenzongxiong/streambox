#ifndef CREATE_H_
#define CREATE_H_

#include <vector>
#include "Values.h"
#include "core/Transforms.h"
#include "Source/Bounded.h"

using namespace std;

//////////////////////////////////////////////

/**
 * A {@code PTransform} that creates a {@code PCollection} from a set of in-memory objects.
 */
template<class T>
class Values: public PTransform {

public:
	PCollection *apply(PValue *input);

private:
	/** The elements of the resulting PCollection. */
	vector<T> _elems;

	// we actually do copy ...
	Values(const vector<T>& elems);

public:

};

//////////////////////////////////////////////

//////////////////////////////////////////////

template <class T> class CreateSource;

template <class T>
class ElementReader: public OffsetBasedReader<T> {
private:
	int _index;
	T _next;

public:
	ElementReader(CreateSource<T>* source)
		: OffsetBasedReader<T>(source), _index(-1) { }

	~ElementReader() { }

	T getCurrent() { return _next; }

	void close() { }

protected:
	long getCurrentOffset() {
		return _index;
	}

	bool startImpl() {
		return advanceImpl();
	}

	CreateSource<T>* getCurrentSource() {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);
		return (CreateSource<T> *) OffsetBasedReader<T>::getCurrentSource();
	}

	bool advanceImpl() {
		CreateSource<T>* source = getCurrentSource();
		_index++;
		if (_index >= (int)(source->_allElements.size())) {
			// XXX we don't use @Optional
			return false;
		}
		_next = source->_allElements[_index];
		return true;
	}
};

//////////////////////////////////////////////


// xzl -- the source for the "Create" transform
// note that we don't store raw bytes nor deal with coder
template <class T>
class CreateSource: public OffsetBasedSource<T> {

public:
  vector<T> _allElements;  // xzl: we own the copy; it is not a pointer.

private:
	long _totalSize;   // xzl: in bytes?

	/**
	 * Create a new source with the specified bytes. The new source owns the input element bytes,
	 * which must not be modified after this constructor is called.
	 *
	 * xzl: note we create a private copy of all elements.
	 */
	CreateSource(const vector<T>& allElements, long totalSize) :
	  OffsetBasedSource<T>::OffsetBasedSource(0, totalSize, 1),
	  _allElements(allElements), _totalSize(totalSize)     { }

	long getEstimatedSizeBytes(PipelineOptions* options) {
		return _totalSize;
	}

	bool producesSortedKeys(PipelineOptions* options) {
		return false;
	}

public:

	static CreateSource* fromIterable(const vector<T>& elements) {
		return new CreateSource(elements, sizeof(T) * elements.size());
	}

	BoundedReader<T>* createReader(PipelineOptions *opt) {
		return new ElementReader<T>(this);
	}

	void validate() {
	}

	long getMaxEndOffset(PipelineOptions* options) {
		return _allElements.size() * sizeof(T);
	}

	// we create duplicated copies of the elements (expensive??)
	// note that this does not shrink this->_allElements.
	// we probably should do splice
	OffsetBasedSource<T>* createSourceForSubrange(long start, long end) {
		// xzl: copy ctor
		vector<T> primaryElems(_allElements.begin() + start / sizeof(T),
				_allElements.end() + end / sizeof(T));
		long primarySizeEstimate = (long) (_totalSize * primaryElems.size()
				/ (double) _allElements.size());

		// xzl: here we make copy again...
		return new CreateSource<T>(primaryElems,
				primarySizeEstimate);
	}

	// xzl: deal with this??
	long getBytesPerOffset() {
		if (_allElements.size() == 0) {
			return 1L;
		}
		return max(1, _totalSize / _allElements.size());
	}

};

//////////////////////////////////////////////

/**
 * Values implementation
 */
template<class T>
PCollection* Values<T>::apply(PValue *input) {
	CreateSource<T>* source = CreateSource<T>::fromIterable(_elems);
	assert(0); // TODO
//		return input->getPipeline()->apply();
	return nullptr;
}

template<class T>
Values<T>::Values(const vector<T>& elems) :
		_elems(elems) { }

//////////////////////////////////////////////

//////////////////////////////////////////////

/**
 * A {@code PTransform} that creates a {@code PCollection} whose elements have
 * associated timestamps.
 *
 * xzl: XXX need to implement DoFn first.
 */
template <class T>
class TimestampedValues: public PTransform {
public:
	PCollection *apply(PValue *input) {
	  return nullptr;
	}

private:
	list<TimestampedValue<T>>* _elems = NULL;

	TimestampedValues(list<TimestampedValue<T>>* elems) {
		_elems = new list<TimestampedValue<T>>(*elems);
	}
};

#endif /* CREATE_H_ */
