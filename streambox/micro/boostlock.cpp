/*
 * boostlock.cpp
 *
 *  Created on: Dec 7, 2016
 *      Author: xzl
 *
 * based on https://gist.github.com/exallium/1976371
 *
 * g++-5 boostlock.cpp -lboost_system -lboost_thread -o boostlock.bin
 */

#include <stdio.h>
#include <iostream>
#include <boost/thread.hpp>

using namespace std;
using boost::thread;
using boost::mutex;
using boost::shared_lock;
using boost::shared_mutex;
using boost::upgrade_lock;
using boost::upgrade_to_unique_lock;

mutex a;
shared_mutex b;

const int THREAD_COUNT = 10;

void worker(int i) {
    a.lock();
    cout << "Worker thread " << i << endl;
    a.unlock();

    shared_lock<shared_mutex> lock(b);
    a.lock();
    cout << "Worker thread " << i << endl;
    a.unlock();
}

int main(int argc, char** argv)
{
    thread *threads[THREAD_COUNT];

    upgrade_lock<shared_mutex> lock(b);
    cout << "before upgrade\n";
    if (lock.owns_lock()) {
    	printf("shared lock is owned\n");
    } else
    	printf("shared lock is NOT owned\n");


    upgrade_to_unique_lock<shared_mutex> uniqueLock(lock);

    cout << "after upgrade\n";
    if (lock.owns_lock()) {
    	printf("shared lock is owned\n");
    } else
    	printf("shared lock is NOT owned\n");

    if (uniqueLock.owns_lock()) {
    	printf("uniqueLock lock is owned\n");
    } else
    	printf("uniqueLock lock is NOT owned\n");

//    // not working
//    upgrade_to_unique_lock<shared_mutex> mylock;
//    mylock = uniqueLock; // try to assign

    // Creation
    for(int i = 0; i < THREAD_COUNT; i++) {
        threads[i] = new thread(worker, i);
    }

    cin.get();
    cout << "Unlocking..." << endl;
    uniqueLock.~upgrade_to_unique_lock();

    // Cleanup
    for(int i = 0; i < THREAD_COUNT; i++) {
        threads[i]->join();
        delete threads[i];
    }

    return 0;
}
