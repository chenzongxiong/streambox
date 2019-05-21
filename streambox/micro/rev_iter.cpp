/*
 * rev_iter.cpp
 *
 *  Created on: Nov 8, 2016
 *      Author: xzl
 *

 g++-5 -std=c++11 rev_iter.cpp -g -O2 -D_GLIBCXX_DEBUG -o rev_iter.bin

 */


/* test the reverse iterator */

#include <iostream>
#include <list>

using namespace std;

void buggy()
{

	list<int> l;

	/* populating... */
	for (int i = 0; i < 10; i ++) {
		l.push_back(i);
	}

	cout << "reverse iteration w/ dump ... \n";
	for (auto it = l.rbegin(); it != l.rend(); it ++) {
		cout << *it << " ";
	}
	cout << endl;

	cout << "reverse iteration /w erase ... \n";
	{
		auto it = l.rbegin();
		while (it != l.rend()) {
			if (*it % 2 == 0) {
				/* bug: @it becomes invalid if we don't reassign it value */
				std::advance(it, 1);
				l.erase(it.base());
			} else
				it ++;
		}
	}

}

int main()
{
	list<int> l;

	/* populating... */
	for (int i = 0; i < 10; i ++) {
		l.push_back(i);
	}

	cout << "forward iteration ... \n";
	for (auto it = l.begin(); it != l.end(); it ++) {
		cout << *it << " ";
	}
	cout << endl;


	cout << "reverse iteration w/ dump ... \n";
	for (auto it = l.rbegin(); it != l.rend(); it ++) {
		cout << *it << " ";
	}
	cout << endl;

	cout << "reverse iteration /w erase ... \n";
	{
		auto it = l.rbegin();
		while (it != l.rend()) {
			if (*it % 2 == 0) {
				it ++;
				it = decltype(l)::reverse_iterator(l.erase(it.base()));
			} else
				it ++;
		}
	}

	cout << "reverse iteration w/ dump ... \n";
	for (auto it = l.rbegin(); it != l.rend(); it ++) {
		cout << *it << " ";
	}
	cout << endl;

	buggy();
}
