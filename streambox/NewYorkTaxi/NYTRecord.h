#ifndef NYTRECORD_H
#define NYTRECORD_H

#include <iostream>
#include <cstdint>
#include <cstring>
#include "Values.h"

struct __attribute__((packed)) NYTRecord {

	char medallion[30];         /* 30 */
	char hack_license[30];      /* 30 */
	char vendorID[4];           /* 4 */
	uint64_t rate_code_id;      /* 8 */
	time_t pickup_ts;           /* 8 */
	time_t dropoff_ts;          /* 8 */
	uint64_t passenger_count;   /* 8 */
	uint64_t trip_time_in_secs; /* 8 */
	double trip_distance;       /* 8 */
	double dropoff_long;        /* 8 */
	double dropoff_lang;        /* 8 */
	double pickup_long;         /* 8 */
	double pickup_lang;         /* 8 */
    uint8_t store_and_fwd_flag; /* 1 */

    NYTRecord(string_range rawRecord) {
        memcpy(&medallion, rawRecord.data, 30);
        memcpy(&hack_license, rawRecord.data + 30, 30);
        memcpy(&vendorID, rawRecord.data + 60, 4);
        memcpy(&rate_code_id, rawRecord.data + 64, 8);
        memcpy(&pickup_ts, rawRecord.data + 72, 8);
        memcpy(&dropoff_ts, rawRecord.data + 80, 8);
        memcpy(&passenger_count, rawRecord.data + 88, 8);
        memcpy(&trip_time_in_secs, rawRecord.data + 96, 8);
        memcpy(&trip_distance, rawRecord.data + 104, 8);
        memcpy(&pickup_long, rawRecord.data + 112, 8);
        memcpy(&pickup_lang, rawRecord.data + 120, 8);
        memcpy(&dropoff_long, rawRecord.data + 128, 8);
        memcpy(&dropoff_lang, rawRecord.data + 136, 8);
        memcpy(&store_and_fwd_flag, rawRecord.data + 144, 1);
    }

    void print() {
        std::cout << std::hex;
		std::cout << " medallion=0x" << medallion;
		std::cout << " hack_license=0x" << hack_license;
		std::cout << " vendorID=0x" << vendorID;
		std::cout << " rate_code_id=0x" << rate_code_id;
		std::cout << " store_and_fwd_flag=0x" << store_and_fwd_flag;
		std::cout << " pickup_ts=0x" << pickup_ts;
		std::cout << " dropoff_ts=0x" << dropoff_ts;
		std::cout << " passenger_count=0x" << passenger_count;
		std::cout << " trip_time_in_secs=0x" << trip_time_in_secs;
		std::cout << " trip_distance=" << trip_distance;
		std::cout << " pickup_long=" << pickup_long;
		std::cout << " pickup_lang=" << pickup_lang;
		std::cout << " dropoff_long=" << dropoff_long;
		std::cout << " dropoff_lang=" << dropoff_lang;
		std::cout << std::endl;
    }
};

#endif
