/*
 * unique_ptr.cpp
 *
 *  Created on: Dec 7, 2016
 *      Author: xzl
 *
 *  Purdue University, 2016
 *
 *  http://www.cplusplus.com/reference/memory/unique_ptr/operator=/
 *
 *  g++-5 unique_ptr.cpp -o unique_ptr.bin -std=c++11
 */

#include <stdio.h>
// unique_ptr::operator= example
#include <iostream>
#include <memory>

void func(std::unique_ptr<int> * pp) {
	*pp = std::move(std::unique_ptr<int>(new int(3)));
//	*pp = std::move(std::unique_ptr<>(new int(3)));
}

int main () {
  std::unique_ptr<int> foo;
  std::unique_ptr<int> bar;
  std::unique_ptr<int> x;

  foo = std::unique_ptr<int>(new int (101));  // rvalue

  bar = std::move(foo);                       // using std::move

  std::cout << "foo: ";
  if (foo) std::cout << *foo << '\n'; else std::cout << "empty\n";

  std::cout << "bar: ";
  if (bar) std::cout << *bar << '\n'; else std::cout << "empty\n";

  std::unique_ptr<int>* p = &bar;
  x = std::move(*p);

  std::cout << "bar: ";
	if (bar) std::cout << *bar << '\n'; else std::cout << "empty\n";

	std::cout << "x: ";
	  if (x) std::cout << *x << '\n'; else std::cout << "empty\n";

	func(&x);

	std::cout << "x: ";
	  if (x) std::cout << *x << '\n'; else std::cout << "empty\n";

  return 0;
}
