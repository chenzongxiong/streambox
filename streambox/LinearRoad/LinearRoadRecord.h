#ifndef LR_RECORD_H
#define LR_RECORD_H

#include <string>
#include <sstream>
#include "Values.h"

struct LinearRoadRecord
{
	int	m_iType;   // Type:
      				   //	. 0: position report
			       	   //	. 2: account balance request
      				   //	. 3: daily expenditure request
			      	   //	. 4: travel time request
	int	m_iTime;   // 0...10799 (second), timestamp position report emitted
	int	m_iVid;	   // 0...MAXINT, vehicle identifier
	int	m_iSpeed;  // 0...100, speed of the vehicle
	int	m_iXway;   // 0...L-1, express way
	int	m_iLane;   // 0...4, lane
	int	m_iDir;    // 0..1, direction
	int	m_iSeg;    // 0...99, segment
	int	m_iPos;    // 0...527999, position of the vehicle
	int	m_iQid;    // query identifier
	int m_iSinit;  // start segment
	int	m_iSend;   // end segment
	int	m_iDow;    // 1..7, day of week
	int	m_iTod;    // 1...1440, minute number in the day
	int	m_iDay;    // 1..69, 1: yesterday, 69: 10 weeks ago

	LinearRoadRecord()	{
		m_iType		= 0;
		m_iTime		= 0;   // 0...10799 (second), timestamp position report emitted
		m_iVid		= 0;	   // 0...MAXINT, vehicle identifier
		m_iSpeed	= 0;  // 0...100, speed of the vehicle
		m_iXway		= 0;	// 0...L-1, express way
		m_iLane		= 0;   // 0...4, lane
		m_iDir		= 0;    // 0..1, direction
		m_iSeg		= 0;    // 0...99, segment
		m_iPos		= 0;    // 0...527999, position of the vehicle
		m_iQid		= 0;    // query identifier
		m_iSinit	= 0;  // start segment
		m_iSend		= 0;   // end segment
		m_iDow		= 0;    // 1..7, day of week
		m_iTod		= 0;    // 1...1440, minute number in the day
		m_iDay		= 0;    // 1..69, 1: yesterday, 69: 10 weeks ago
	}


  LinearRoadRecord(string_range raw_record) {
    char str[128];
    memset(str, '\0', 128);

    strncpy(str, raw_record.data, (size_t)(raw_record.len));
    // zxchen: initialize function is expensive, drop 10 times. why ?
    // simply due to memory movement?
    initialize(str);
    // zxchen: TOO EXPENSIVE, but why ?
    // istringstream ss(raw_record.data);
    // std::vector<std::string> tokens;
    // std::string token;
    // while (getline(ss, token, ',')) {
    //     tokens.push_back(token);

    // }
    // if (tokens.size() != 15) {
    //     std::cout << "The size: " << tokens.size() << " is invalid" << std::endl;
    //     abort();
    // }
    // // int i = 0;
    // // for (auto it = tokens.begin(); it != tokens.end(); it ++) {
    // //     std::cout << "token[" << i++ << "]: " << *it << "\t";
    // // }
    // // std::cout << endl;
    // m_iType = std::stoi(tokens[0]);
    // m_iTime= std::stoi(tokens[1]);
    // m_iVid = std::stoi(tokens[2]);
    // m_iSpeed = std::stoi(tokens[3]);
    // m_iXway = std::stoi(tokens[4]);
    // m_iLane = std::stoi(tokens[5]);
    // m_iDir = std::stoi(tokens[6]);
    // m_iSeg = std::stoi(tokens[7]);
    // m_iPos = std::stoi(tokens[8]);
    // m_iQid = std::stoi(tokens[9]);
    // m_iSinit = std::stoi(tokens[10]);
    // m_iSend = std::stoi(tokens[11]);
    // m_iDow = std::stoi(tokens[12]);
    // m_iTod = std::stoi(tokens[13]);
    // m_iDay = std::stoi(tokens[14]);
    if (! valid()) {
        // roughtly check
        std::cout << raw_record.data << std::endl;
        std::cout << toString() << std::endl;
    }
  }

	void initialize(char* str) {
        char *token = NULL;

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iType = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iTime = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iVid = atoi(token);	}

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iSpeed = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iXway = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iLane = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iDir = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iSeg = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iPos = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iQid = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iSinit = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iSend = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iDow = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iTod = atoi(token); }

        token = strtok_r( str, ",", &str);
		if (token != NULL) {m_iDay = atoi(token); }


	}

  std::string toString() {
    std::ostringstream oss;
    oss << "Type: " << m_iType << "\tTime: " << m_iTime << "\tVid: " << m_iVid
        << "\tSpeed: " << m_iSpeed << "\tXway: " << m_iXway << "\tLane: " << m_iLane
        << "\tDir: " << m_iDir << "\tSeg: " << m_iSeg << "\tPos: " << m_iPos
        << "\tQid: " << m_iQid << "\tSinit: " << m_iSinit << "\tSend: " << m_iSend
        << "\tDow: " << m_iDow << "\tTod: " << m_iTod << "\tDay: " << m_iDay;
    return oss.str();
  }
    bool valid() {
        return (m_iType >= -1) && (m_iTime >= -1) && (m_iVid >= -1) && (m_iSpeed >= -1) &&
            (m_iXway >= -1) && (m_iLane >= -1) && (m_iDir >= -1) && (m_iSeg >= -1) &&
            (m_iQid >= -1) && (m_iSinit >= -1) && (m_iSend >= -1) && (m_iDow >= -1) &&
            (m_iTod >= -1) && (m_iDay >= -1);
    }
};

#endif /* LR_RECORD_H */
