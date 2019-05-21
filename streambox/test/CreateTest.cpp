#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <typeinfo>

#ifdef USE_CILK
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#endif

#include "Create.h"
#include "Read.h"

//#define xstr(s) str(s)
//#define str(s) #s
#define xstr(s) #s //hym: avoid conflict with .str()

/* expect: Pipeline* p = Pipeline::create(NULL) */
#define source_transform(T) \
	  auto T##_output = dynamic_cast<PCollection *>(p->apply1(&T)); \
	  T##_output->_name = xstr(T##_out)
	  //T##_output->_name = str(T##_out)

/* only works for apply1() */
#define connect_transform(T1, T2) \
    auto T2##_output = dynamic_cast<PCollection *>(T1##_output->apply1(&T2)); \
	T2##_output -> _name = xstr(T2##_out)
	//T2##_output -> _name = str(T2##_out)


void testRead()
{
	vector<long> v ({1,2,3,4,5,6,7,8,9}); // 10x longs
	CreateSource<long>* source = \
		CreateSource<long>::fromIterable(v);

	// create a new transform
	Bounded<long>* b = Read<long>::from(source);
	cout << b->getKindString();
}

void testCreateSource()
{
	vector<long> v ({1,2,3,4,5,6,7,8,9}); // 10x longs
	CreateSource<long>* source = CreateSource<long>::fromIterable(v);
	auto reader = new ElementReader<long>(source);
	reader->start();
	for (int i = 0; i < 20; i++) {
		printf("value is %ld, timestamp %s\n",
				reader->getCurrent(),
				to_simple_string(reader->getCurrentTimestamp()).c_str());

		if (!reader->advance())
			break;
	}
}

#include "core/Pipeline.h"
void testPipeline()
{
	vector<long> v ({1,2,3,4,5,6,7,8,9}); // 10x longs
	CreateSource<long>* source = \
			CreateSource<long>::fromIterable(v);

	// create a new Bounded transform
	Bounded<long>* b = Read<long>::from(source);
	cout << b->getKindString();

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	// this will finally go to Bounded::apply() which returns a PTransform*
	PValue *o = p->apply(b);

	cout << "output has type: " << typeid(*o).name() << endl;

	auto c = dynamic_cast<PCollection *>(o);

	// try a transform
	PTransformNop<long> nop;
	c->apply(&nop);

	// try a visitor (the basis of evaluation)
	PipelineVisitorNop visitor;
	p->traverseTopologically(&visitor);
}

//#include "TransformEvaluator.h"
#include "TransformEvaluate.h"
void testParDo()
{
	vector<long> v (100, 42); //
	CreateSource<long>* source = \
				     CreateSource<long>::fromIterable(v);

	// create a new Bounded transform
	Bounded<long>* b = Read<long>::from(source);
	cout << b->getKindString();

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	// this will finally go to Bounded::apply() which returns a PTransform*
	PValue *o = p->apply(b);

	auto c = dynamic_cast<PCollection *>(o);

	// NOP transform
	PTransformNop<long> nop;
	auto c1 = c->apply(&nop);

	// ParDo transform
	DoFnPlus1 dofn;
	ParDo<long, long> pardo("pardo", &dofn);
	dynamic_cast<PCollection *>(c1)->apply(&pardo);

	// Eval the pipeline
	EvaluationBundleContext eval;
	eval.run(p); // this should traverse the pipeline

	printf("eval pipeline done\n");
}

// the "bundle" version
void testParDo1()
{
	vector<long> v (100, 42); // 100x longs
	CreateSource<long>* source = \
				     CreateSource<long>::fromIterable(v);

	// create a new Bounded transform
	Bounded<long>* b = Read<long>::from(source, "bounded");

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PCollection *c = dynamic_cast<PCollection *>(p->apply1(b));

	// NOP transform
	PTransformNop<long> nop;
	auto c1 = dynamic_cast<PCollection *>(c->apply1(&nop));

	DoFnPrinter<long> dofn;
	ParDo<long, long> pardo("pardo-printer", &dofn);
	c1->apply1(&pardo);

	// Eval the pipeline
	EvaluationBundleContext eval;
	eval.runSimple(p);

	printf("eval pipeline done\n");
}

#include <tuple>
void testParDoTuple()
{
	typedef tuple<long, float> Tuple;

	// element type: (long, float)
	vector<Tuple> v ({
			std::make_tuple(1, 10.00),
			std::make_tuple(2, 20.00),
			std::make_tuple(3, 30.00),
			std::make_tuple(4, 40.00),
			std::make_tuple(5, 50.00),
			});

	CreateSource<Tuple>* source = \
				      CreateSource<Tuple>::fromIterable(v);

	// create a new Bounded transform
	Bounded<Tuple>* b = Read<Tuple>::from(source);
	cout << b->getKindString();

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	// this will finally go to Bounded::apply() which returns a PTransform*
	PValue *o = p->apply(b);

	auto c = dynamic_cast<PCollection *>(o);
	assert(c);

	// NOP transform
	PTransformNop<Tuple> nop;
	auto c1 = c->apply(&nop);

	// ParDo transform
	DoFnPlus1Tuple2<long, float> dofn;
	ParDo<Tuple, Tuple> pardo("pardo", &dofn);
	dynamic_cast<PCollection *>(c1)->apply(&pardo);

	// Eval the pipeline
	EvaluationBundleContext eval;
	eval.run(p); // this should traverse the pipeline

	printf("eval pipeline done\n");
}

void testParDoCol()
{
	typedef tuple<long, float> Tuple;

	// element type: (long, float)
#if 0
	vector<Tuple> v ({
			std::make_tuple(1, 10.00),
			std::make_tuple(2, 20.00),
			std::make_tuple(3, 30.00),
			std::make_tuple(4, 40.00),
			std::make_tuple(5, 50.00),
			});
#endif

	vector<Tuple> v (10, std::make_tuple(1, 10.00));

	CreateSource<Tuple>* source = \
				      CreateSource<Tuple>::fromIterable(v);

	// create a new Bounded transform
	BoundedCol2<long, float>* b = ReadCol2<long, float>::from(source, "bounded");

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	// this will finally go to Bounded::apply() which returns a PTransform*
	PValue *o = p->apply(b);

	auto c = dynamic_cast<PCollection *>(o);
	assert(c);

	// Note: cannot use scope here, as objects will be destroyed upon pipeline
	// execution.

#if 1
	// ParUpdateCol2 transform
	DoFnSimplePlus1<long> fn1;
	DoFnSimplePlus1<float> fn2;

	auto parupdate = new ParUpdateCol2<long, float> ("parupdatecol", &fn1, &fn2);
	auto c2 = c->apply(parupdate);
	assert(c2);
#endif

#if 1
	DoFnSimplePrinter<long> fn11;
	DoFnSimplePrinter<float> fn22;

	ParUpdateCol2<long, float> printer ("printcol", &fn11, &fn22);
	dynamic_cast<PCollection *>(c2)->apply(&printer);
#endif

	// Eval the pipeline
	EvaluationBundleContext eval;
	eval.run(p); // this should traverse the pipeline

	printf("eval pipeline done\n");
}

void testUnbounded() {

	// create a new Bounded transform
	UnboundedPseudo<long> b ("unbounded", 100);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	// this will finally go to Bounded::apply() which returns a PTransform*
	PValue *o = p->apply(&b);

	auto c1 = dynamic_cast<PCollection *>(o);

	DoFnPrinter<long> dofn;
	ParDo<long, long> pardo("pardo-printer", &dofn);
	c1->apply(&pardo);

	// Eval the pipeline
	EvaluationBundleContext eval;
	eval.run(p);

	printf("eval pipeline done\n");
}

void testFixedWindowInto() {

	// create a new Bounded transform
	UnboundedPseudo<long> b ("unbounded", 50);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	// this will finally go to Bounded::apply() which returns a PTransform*
	PValue *o = p->apply(&b);

	auto c1 = dynamic_cast<PCollection *>(o);
	c1->_name = "src_out";

	FixedWindowInto<long> fwi ("window", seconds(3));
	auto c2 = dynamic_cast<PCollection *>(c1->apply(&fwi));
	c2->_name = "win_out";

	WinSum<long> wsum("sum");
	auto c3 = dynamic_cast<PCollection *>(c2->apply(&wsum));
	c3->_name = "sum_out";

	DoFnPrinter<long> dofn;
	ParDo<long, long> pardo("pardo-printer", &dofn);
	auto c4 = dynamic_cast<PCollection *>(c3->apply(&pardo));
	c4->_name = "printer_out";

	// Eval the pipeline
	EvaluationBundleContext eval;
	eval.run(p);

	printf("eval pipeline done\n");
}

void testSelect() {

	// create a new Bounded transform
	UnboundedPseudo<long> b ("unbounded", 50);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PValue *o = p->apply(&b);

	auto c1 = dynamic_cast<PCollection *>(o);
	c1->_name = "src_out";

	SimpleSelect<long> sel ("select");
	auto c11 = dynamic_cast<PCollection *>(c1->apply(&sel));
	c11->_name = "select_out";

	FixedWindowInto<long> fwi ("window", seconds(3));
	auto c2 = dynamic_cast<PCollection *>(c11->apply(&fwi));
	c2->_name = "win_out";

	WinSum<long> wsum("sum");
	auto c3 = dynamic_cast<PCollection *>(c2->apply(&wsum));
	c3->_name = "sum_out";

	DoFnPrinter<long> dofn;
	ParDo<long, long> pardo("pardo-printer", &dofn);
	auto c4 = dynamic_cast<PCollection *>(c3->apply(&pardo));
	c4->_name = "printer_out";

	// Eval the pipeline
	EvaluationBundleContext eval;
	eval.run(p);

	printf("eval pipeline done\n");
}

/*
 * Unbounded
 *     |  RecordBitmapBundle<long>
 *     V
 * SimpleMapper
 *     |  RecordBitmapBundle<long,long>
 *     V
 *   WinGBK
 *     |  WindowsKeyedBundle<KVPair>
 *     V
 *  WindowKeyedReducer
 *     | WindowsKeyedBundle<KVPair>
 *     V
 */
void testWinGBK()
{
	// create a new Bounded transform
	UnboundedPseudo<long> b ("[unbounded]", 50);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PValue *o = p->apply(&b);

	auto c1 = dynamic_cast<PCollection *>(o);
	c1->_name = "src_out";

	SimpleMapper<long, pair<long, long>> mapper ("[mapper]");
	//  auto c11 = dynamic_cast<PCollection<pair<long,long>> *>(c1->apply(&mapper));
	auto c11 = dynamic_cast<PCollection *>(c1->apply(&mapper));
	c11->_name = "mapper_out";

	WinGBK<pair<long, long>> wgbk ("[wingbk]", seconds(1));
	auto c2 = dynamic_cast<PCollection *>(c11->apply(&wgbk));
	c2->_name = "wingbk_out";

	WindowKeyedReducer<pair<long,long>> reducer("[reducer]");
	auto c3 = dynamic_cast<PCollection *>(c2->apply(&reducer));
	c3->_name = "reducer";

#if 0
	DoFnPrinter<pair<long,long>> dofn;
	ParDo<pair<long,long>, pair<long,long>> pardo("pardo-printer", &dofn);
	auto c4 = dynamic_cast<PCollection *>(c2->apply(&pardo));
	c4->_name = "printer_out";
#endif

	// Eval the pipeline
	EvaluationBundleContext eval;
	eval.run(p);

	printf("eval pipeline done\n");
}

#if 0
void testJoin()
{
	// create a new Bounded transform
	UnboundedPseudo<long> unbound ("[unbounded]", 50);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	vector<PCollection*> o2 = p->apply2(&unbound);
	o2[0]->_name = "src_out0";
	o2[1]->_name = "src_out1";

	SimpleMapper<long, pair<long, long>> mapper0 ("[mapper0]");
	auto c11 = dynamic_cast<PCollection *>(o2[0]->apply1(&mapper0));
	c11->_name = "mapper0_out";

	SimpleMapper<long, pair<long, long>> mapper1 ("[mapper1]");
	auto c12 = dynamic_cast<PCollection *>(o2[1]->apply1(&mapper1));
	c12->_name = "mapper1_out";

#if 1
	Join<pair<long,long>> join("join", seconds(3));
	auto c2 = dynamic_cast<PCollection *>(c11->apply1(&join));
	auto c3 = dynamic_cast<PCollection *>(c12->apply1(&join));
	assert(c2 == c3); // same output.
#endif

	EvaluationBundleContext eval;
	eval.runSimple(p);
}
#endif

/*
 *                    +--- SimpleMapper -- (KVPair)--\
 * Unbound -(long)->  +                             Join --> (KVPair)
 *                    +--- SimpleMapper -- (KVPair)--/
 *
 */
void testJoin()
{
	// create a new Bounded transform
	UnboundedInMem<long> unbound ("[unbounded]", 50);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	vector<PCollection*> o2 = p->apply2(&unbound);
	o2[0]->_name = "src_out0";
	o2[1]->_name = "src_out1";

	SimpleMapper<long, pair<long, long>> mapper0 ("[mapper0]");
	auto c11 = dynamic_cast<PCollection *>(o2[0]->apply1(&mapper0));
	c11->_name = "mapper0_out";

	SimpleMapper<long, pair<long, long>> mapper1 ("[mapper1]");
	auto c12 = dynamic_cast<PCollection *>(o2[1]->apply1(&mapper1));
	c12->_name = "mapper1_out";

#if 1
	Join<pair<long,long>> join("join", seconds(3));
	auto c2 = dynamic_cast<PCollection *>(c11->apply1(&join));
	auto c3 = dynamic_cast<PCollection *>(c12->apply1(&join));
	assert(c2 == c3); // same output.
#endif

	RecordBitmapBundleSink<pair<long, vector<long>>> sink("[sink]");
	auto c4 = dynamic_cast<PCollection *>(c3->apply1(&sink));
	assert(c4); //just for fixing warning

	EvaluationBundleContext eval;
	eval.runSimple(p);
}

void testSimplePipeline()
{
	// create a new Bounded transform
	UnboundedPseudo<long> unbound ("[unbounded]", 50);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	vector<PCollection*> o2 = p->apply2(&unbound);
	o2[0]->_name = "src_out0";
	o2[1]->_name = "src_out1";

	SimpleMapper<long, pair<long, long>> mapper0 ("[mapper0]");
	auto c11 = dynamic_cast<PCollection *>(o2[0]->apply1(&mapper0));
	c11->_name = "mapper0_out";

	SimpleMapper<long, pair<long, long>> mapper1 ("[mapper1]");
	auto c12 = dynamic_cast<PCollection *>(o2[1]->apply1(&mapper1));
	c12->_name = "mapper1_out";

	EvaluationBundleContext eval;
	eval.runSimple(p);
}

	template<class T>
auto operator<<(std::ostream& os, const T& t) -> decltype(t.print(os), os)
{
	t.print(os);
	return os;
}

void testSessionWindowInto()
{
	using date = boost::gregorian::date;

	{
		W("test1: inorder record insertion");
		// each record has a 3-sec window
		SessionWindowInto<long> trans("swi", seconds(3));

		vector<Record<long>> recs;
		ptime current (date(2002,Jan,10));

		for (int i = 0; i < 100; i++) {
			if (!(i % 10))
				current += seconds(10); // leap
			else
				current += seconds(1);

			recs.emplace_back(i, current);
		}

		for (auto && it = recs.begin(); it != recs.end(); it++) {
			trans.try_add_record(&(*it), nullptr);
		}

		trans.dump();
	}

	{
		W("test2: window merge");

		// each record has a 3-sec window
		SessionWindowInto<long> trans("swi", seconds(3));
		ptime start (date(2002,Jan,10));

		vector<Record<long>> recs (
				{
				Record<long>(1, start + seconds(0)),
				Record<long>(2, start + seconds(4)),
				Record<long>(3, start + seconds(8))
				}
				);

		for (auto && it = recs.begin(); it != recs.end(); it++) {
			trans.try_add_record(&(*it), nullptr);
		}
		trans.dump();


		W("test2: after merging...");
		vector<Record<long>> recs2 (
				{
				Record<long>(4, start + seconds(2)),
				Record<long>(5, start + seconds(6)),
				}
				);

		for (auto && it = recs2.begin(); it != recs2.end(); it++) {
			trans.try_add_record(&(*it), nullptr);
		}
		trans.dump();
	}


	{
		W("test3: reverse order record insertion");
		// each record has a 3-sec window
		SessionWindowInto<long> trans("swi", seconds(3));

		vector<Record<long>> recs;
		ptime current (date(2002,Jan,10));

		for (int i = 0; i < 100; i++) {
			if (!(i % 10))
				current -= seconds(10); // leap
			else
				current -= seconds(1);

			recs.emplace_back(i, current);
		}

		for (auto && it = recs.begin(); it != recs.end(); it++) {
			trans.try_add_record(&(*it), nullptr);
		}

		trans.dump();
	}

	{
		W("test4: window merge (slow path)");

		// each record has a 3-sec window
		SessionWindowInto<long> trans("swi", seconds(3));
		ptime start (date(2002,Jan,10));

		vector<Record<long>> recs (
				{
				Record<long>(1, start + seconds(10)),
				Record<long>(2, start + seconds(14)),
				Record<long>(3, start + seconds(18))
				}
				);

		for (auto && it = recs.begin(); it != recs.end(); it++) {
			trans.try_add_record(&(*it), nullptr);
		}
		trans.dump();


		W("test4: after merging...");
		vector<Record<long>> recs2 (
				{
				Record<long>(4, start + seconds(9)),
				Record<long>(5, start + seconds(13)),
				}
				);

		for (auto && it = recs2.begin(); it != recs2.end(); it++) {
			trans.try_add_record(&(*it), nullptr);
		}
		trans.dump();
	}

	{
		W("test5: handle bundle input");

		// each record has a 3-sec window
		SessionWindowInto<long> trans("swi", seconds(3));
		auto bundle = make_shared<RecordBundleDebug<long>>();

		ptime start (date(2002,Jan,10));

		vector<Record<long>> recs (
				{
				Record<long>(1, start + seconds(10)),
				Record<long>(2, start + seconds(14)),
				Record<long>(3, start + seconds(18))
				}
				);

		for (auto && it = recs.begin(); it != recs.end(); it++) {
			bundle->add_record(*it);
		}

		trans.MergeRecordBundle(bundle);
		trans.dump();

		W("test5: retrieve windows");
		auto winmap = trans.RetrieveWindows(true, start + seconds(15));
		cout << trans;

		for (auto &  kv : winmap) {
			cout << (*kv.second);
		}

		cout << "---- advance watermark and close/purge the remaining windows ---- \n";

		winmap = trans.RetrieveWindows(true, start + seconds(100));
		cout << trans;

		for (auto &  kv : winmap) {
			cout << (*kv.second);
		}

		cout << " ---- all done ---- \n";
	}

	{
		W("test6: merge window sets");

		// each record has a 3-sec window
		SessionWindowInto<long> trans1("swi", seconds(3));
		SessionWindowInto<long> trans2("swi", seconds(3));

		auto bundle1 = make_shared<RecordBundleDebug<long>>();
		auto bundle2 = make_shared<RecordBundleDebug<long>>();

		ptime start (date(2002,Jan,10));

		vector<Record<long>> recs (
				{
				Record<long>(1, start + seconds(10)),
				Record<long>(2, start + seconds(14)),
				Record<long>(3, start + seconds(18))
				}
				);

		vector<Record<long>> recs2 (
				{
				/* out of any bundle1's session window */
				Record<long>(50, start + seconds(1)),
				/* same as bundle1's session windows */
				Record<long>(100, start + seconds(10)),
				Record<long>(200, start + seconds(14)),
				Record<long>(300, start + seconds(18)),
				/* partially overlap with bundle1's windows */
				Record<long>(350, start + seconds(19)),
				/* out of any bundle1's session window */
				Record<long>(400, start + seconds(50))
				}
				);

		for (auto && it = recs.begin(); it != recs.end(); it++) {
			bundle1->add_record(*it);
		}

		for (auto && it = recs2.begin(); it != recs2.end(); it++) {
			bundle2->add_record(*it);
		}

		trans1.MergeRecordBundle(bundle1);
		trans2.MergeRecordBundle(bundle2);
		trans1.dump();

		trans1.MergeWindowsSet(trans2.windows);

		cout << "after merge window sets....\n";
		cout << trans1;

		auto winmap = trans1.RetrieveWindows(true, start + seconds(500));
		cout << trans1;
		cout << "emitting windowed records ...\n";
		for (auto &  kv : winmap) {
			cout << (*kv.second);
		}

		cout << " ---- all done ---- \n";
	}

	{
		W("test7: pipeline");
		// create a new Bounded transform
		UnboundedPseudo<long> unbound ("[unbounded]",
				50, /* bundle interval */
				5000  /* session interval */
				);

		// create a new pipeline
		Pipeline* p = Pipeline::create(NULL);

		PCollection *o2 = dynamic_cast<PCollection *>(p->apply1(&unbound));
		o2->_name = "src_out0";

		SessionWindowInto<long> swi ("swi", seconds(1));
		auto c11 = dynamic_cast<PCollection *>(o2->apply1(&swi));
		c11->_name = "swi_out";

		EvaluationBundleContext eval;
		eval.runSimple(p);
	}

}

/*
 *      UnboundedInMem
 *            | RecordBitmapBundle<string_range>
 *            V
 *       WordCountMapper
 *            | RecordBitmapBundle<KVPair<string,long>>
 *            V
 *          WinGBK
 *            | WindowsKeyedBundle<KVPair<string,long>>
 *            V
 *      WindowKeyedReducer (stateful)
 *            | WindowsKeyedBundle<KVPair<string,long>>
 *            V
 *
 */
void testWordCount()
{
	//  UnboundedInMem<string_range> unbound("unbounded-inmem", "/ssd/word_100MB.txt");
	UnboundedInMem<string_range> unbound("unbounded-inmem",
			//      "/home/xzl/Dropbox/private/source/creek/test/Debug/word_10MB.txt",
			//      "/home/xzl/Dropbox/private/source/creek/test/Debug/big.txt",
			"/home/miaohy/Project/creek/test/Debug/word_100MB.txt",       
			100); // 1000 seems fine

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PCollection *o2 = dynamic_cast<PCollection *>(p->apply1(&unbound));
	o2->_name = "src_out";

	WordCountMapper<> mapper ("[wc-mapper]");
	auto c11 = dynamic_cast<PCollection *>(o2->apply1(&mapper));
	c11->_name = "mapper_out";

	WinGBK<pair<string, long>> wgbk ("[wingbk]", seconds(1));
	auto c2 = dynamic_cast<PCollection *>(c11->apply1(&wgbk));
	c2->_name = "wingbk_out";

	WindowKeyedReducer<pair<string,long>> reducer("[reducer]");
	auto c3 = dynamic_cast<PCollection *>(c2->apply1(&reducer));
	c3->_name = "reducer";

	EvaluationBundleContext eval;
	eval.runSimple(p);
}

/*
 *      UnboundedInMem
 *           |  RecordBitmapBundle<string_range>
 *           V
 *        GrepMapper
 *           |  RecordBitmapBundle<string>
 *           V
 *
 *  NB: this needs more work. Consider using testWnidowedGrep()
 */
void testGrep()
{
	//  UnboundedInMem<string_range> unbound("unbounded-inmem", "/ssd/word_100MB.txt");
	UnboundedInMem<string_range> unbound("unbounded-inmem",
			//      "/home/xzl/Dropbox/private/source/creek/test/Debug/word_10MB.txt",
			//      "/home/xzl/Dropbox/private/source/creek/test/Debug/big.txt",
			"/home/miaohy/Project/creek/test/Debug/word_100MB.txt", 
			100); // 1000 seems fine

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PCollection *o2 = dynamic_cast<PCollection *>(p->apply1(&unbound));
	o2->_name = "src_out";

	GrepMapper<> mapper ("HELLO", "[grep-mapper]");
	auto c11 = dynamic_cast<PCollection *>(o2->apply1(&mapper));
	c11->_name = "mapper_out";

	/* XXX windowing XXX */

	EvaluationBundleContext eval;
	eval.runSimple(p);
}

/*
 *      UnboundedInMem
 *           |  RecordBitmapBundle<string_range>
 *           V
 *      FixedWindowInto
 *           |  WindowsBundle<string_range>
 *           V
 *     WindowedGrepMapper
 *           |  WindowsBundle<string>
 *           V
 *       WindowedSum (stateful)
 *           |  WindowsBundle<vector<string>>
 *           V
 *          Sink
 *
 */
void testWindowedGrep()
{
	UnboundedInMem<string_range> unbound("unbounded-inmem",
			//"/home/xzl/Dropbox/private/source/creek/test/Debug/small.txt",
			"/home/miaohy/Project/creek/test/Debug/word_100MB.txt",
			100);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
	unbound_output->_name = "unbound_output";

	FixedWindowInto<string_range> fwi ("window", seconds(1));
	//		auto c2 = dynamic_cast<PCollection *>(o2->apply1(&fwi));
	//		c2->_name = "win_out";
	connect_transform(unbound, fwi);

	/* four digits as year */
	//	  WindowedGrepMapper<> grep (R"(\d\d\d\d)", "mapper");

	/* url */
	WindowedGrepMapper<> grep (
			R"([-a-zA-Z0-9@:%._\+~#=]{2,256}\.[a-z]{2,6}\b([-a-zA-Z0-9@:%_\+.~#?&//=]*))",
			"mapper");

	connect_transform(fwi, grep);

	WinSum<string, vector<string>> agg ("agg");
	//		auto c3 = dynamic_cast<PCollection *>(c2->apply1(&agg));
	//		c3->_name = "agg_out";
	connect_transform(grep, agg);

	RecordBundleSink<vector<string>> sink("sink");
	connect_transform(agg, sink);

	EvaluationBundleContext eval;
	eval.runSimple(p);
}

using test_func_ptr_t = void (*) ();

struct test_func_t {
	char const * name;
	test_func_ptr_t func;
};

#define declare_test(fname) { \
	.name = #fname, \
	.func = fname \
},

struct test_func_t test_funcs[] = {
	declare_test(testCreateSource)
		declare_test(testRead)
		declare_test(testPipeline)
		declare_test(testParDo)
		declare_test(testParDo1)
		declare_test(testParDoTuple)
		declare_test(testParDoCol)
		declare_test(testUnbounded)
		declare_test(testFixedWindowInto)
		declare_test(testSelect)
		declare_test(testWinGBK)
		declare_test(testSimplePipeline)
		declare_test(testJoin)
		declare_test(testSessionWindowInto)
		declare_test(testWordCount)
		declare_test(testGrep)
		declare_test(testWindowedGrep)
		/* add more tests here */
};

int main(int argc, char **argv)
{
#if 0
#ifdef USE_CILK
	if (0!= __cilkrts_set_param("nworkers","16"))
	{
		printf("Failed to set worker count\n");
		return 1;
	}
#endif
#endif

	int num_tests = sizeof(test_funcs) / sizeof(test_func_t);

	if (argc < 2) {
		printf("choose tests --- (run %s ${testid}) \n", argv[0]);
		for (int i = 0; i < num_tests; i++) {
			auto && f = test_funcs[i];
			printf("%d: %s\n", i, f.name);
		}
	} else {
		int c = atoi(argv[1]);
		if (c == -1)
			c = num_tests - 1;
		printf("run test %d %s ...\n", c, test_funcs[c].name);
		test_funcs[c].func();
	}
}


