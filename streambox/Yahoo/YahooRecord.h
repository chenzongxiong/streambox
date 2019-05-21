#ifndef YAHOORECORD_H
#define YAHOORECORD_H

#include <cstdint>
#include <cstring>
#include "Values.h"

struct __attribute__((packed)) YahooRecord {
  uint8_t user_id[16];
  uint8_t page_id[16];
  uint8_t campaign_id[16];
  char ad_type[9];
  char event_type[9];
  int64_t current_ms;
  uint32_t ip;

	YahooRecord(string_range rawRecord) {
        // TOO EXPENSIVE.
		memcpy(&user_id, rawRecord.data, 16);
		memcpy(&page_id, rawRecord.data + 16, 16);
		memcpy(&campaign_id, rawRecord.data + 32, 16);
		memcpy(&ad_type, rawRecord.data + 48, 9);
		memcpy(&event_type, rawRecord.data + 57, 9);
		memcpy(&current_ms, rawRecord.data + 66, 8);
		memcpy(&ip, rawRecord.data + 74, 4);
    }

    void print() {
        std::cout << "user_id: " << user_id << std::endl;
        std::cout << "page_id: " << page_id << std::endl;
        // std::cout << "user_id: " << std::bitset<128>(user_id) << std::endl;
        // std::cout << "page_id: " << page_id << std::endl;
        // std::cout << "campaign_id:" << campaign_id << std::endl;
        std::cout << "ad_type: " << ad_type << std::endl;
        std::cout << "current_ms: " << std::bitset<64>(current_ms) << std::endl;
        std::cout << "ip: " << std::bitset<32>(ip) << std::endl;
    }
};

#endif /* YAHOORECORD_H */
