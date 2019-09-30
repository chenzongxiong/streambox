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

    // void setPrice(uint64_t p) {
    //     price = p;
    // }
    NexmarkRecord(uint64_t a, uint64_t b, uint64_t p, uint64_t d) {
        auction = a;
        bidder = b;
        price = p;
        dateTime = d;
    }

}; // size 32 bytes


struct __attribute__((packed)) NexmarkOutputRecord {
    uint64_t auction;
    uint64_t price;
    NexmarkOutputRecord() {}
    NexmarkOutputRecord(NexmarkRecord &rec) {
        auction = rec.auction;
        price = rec.price;
    }
    NexmarkOutputRecord(uint64_t auction, uint64_t price) {
        this->auction = auction;
        this->price = price;
    }
};

#endif /* NEXMARK_RECORD_HPP */
