#ifndef NEXMARK_RECORD_HPP
#define NEXMARK_RECORD_HPP
#include "Values.h"

struct __attribute__((packed)) NexmarkRecord {
    uint64_t auction;
    uint64_t bidder;
    uint64_t price;
    uint64_t dateTime;

    NexmarkRecord() {}

    NexmarkRecord(string_range rawRecord) {
        memcpy(&auction, rawRecord.data, 8);
        memcpy(&bidder, rawRecord.data, 8);
        memcpy(&price, rawRecord.data, 8);
        memcpy(&dateTime, rawRecord.data, 8);
    }

    void print() {
        std::cout << "auction: " << auction
                  << "\nbidder: " << bidder
                  << "\nprice: " << price
                  << "\ndateTime: " << dateTime
                  << std::endl;
    }
}; // size 32 bytes

#endif /* NEXMARK_RECORD_HPP */
