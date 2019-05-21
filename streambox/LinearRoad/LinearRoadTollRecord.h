#ifndef LRB_TOLL_RECORD_H
#define LRB_TOLL_RECORD_H

#include <string>
#include <sstream>
#include <map>

#define SEGMENT_LENGTH 5280
#define SEGMENT_NUM 100


struct __attribute__((packed)) StopMap {
  StopMap() {
    pos = -1;
    count = 0;
  }
  int pos;
  int count;
  std::string toString() {
    std::ostringstream oss;
    oss << "[StopMap]::Pos: " << pos << "\tCount: " << count;
    return oss.str();
  }
};

struct __attribute__((packed)) AvgSpeedMap { /*index is the vecId */

  AvgSpeedMap() {
    value = 0;
    count = 0;
  };
  size_t value;
  size_t count;
};

class AccidentDetection {
public:
  struct __attribute__((packed)) accAtPos {
    int vId1;
    int vId2;
    int time;
  };

  AccidentDetection() {
    deleteCallCnt = 0;
    addCallCnt = 0;
    couldNotBe = 0;
    addFirst = 0;
    addSecond = 0;
    addNew = 0;
    addAlreadyFirst = 0;
    addAlreadySecond = 0;
    removeFirst = 0;
    removeSecond = 0;
    erase = 0;

  }

  void addStop(int vID, int pos, int time)
  {
    addCallCnt++;
    //if a stop is already available
    if(m.count(pos) != 0)
    {
      //One car already at that pos
      if(m[pos].vId1 == vID)
      {
        //do nothing
        addAlreadyFirst++;
      }
      else if(m[pos].vId2 == vID)
      {
        //do nothing
        addAlreadySecond++;

      }
      else if(m[pos].vId2 == 0)
      {
        m[pos].vId2 = vID;
        addSecond++;
      }
      else
      {
        //				cout << "vid1=" << m[pos].vId1 << " vid2=" << m[pos].vId2 << endl;
        couldNotBe++;
        return;
        assert(0);
      }
    }
    else
    {
      //add new entry to pos
      m[pos].vId1 = vID;
      addNew++;
    }

    //update time
    m[pos].time = time;

#ifdef detectEachInsert
    detectAccident();
#endif
  }
  bool isInAccidentMapEntryForPos(int pos)
  {
    return m.count(pos) != 0;
  }

  bool isVedAtPos(int vec, int pos)
  {
    if(m[pos].vId1 == vec)
    {
      return true;
    }
    else if(m[pos].vId2 == vec)
    {
      return true;
    }
    else
      return false;
  }

  void removeEntry(int vecID, int pos)
  {
    deleteCallCnt++;
    if(m[pos].vId1 == vecID)
    {
      m[pos].vId1 = 0;
      removeFirst++;

    }
    else if(m[pos].vId2 == vecID)
    {
      m[pos].vId2 = 0;
      removeSecond++;
    }
    else
    {
      assert(0);
    }

    if(m[pos].vId1 == 0 && m[pos].vId2 == 0)
    {
      m.erase(pos);
      erase++;
      //			cout<<  "erease acc"<< endl;
    }
  }

  void removeStop(int vID, int pos)
  {
    //if a stop is already available
    if(m.count(pos) != 0)
    {
      //One car already at that pos
      if(m[pos].vId1 == vID)
      {
        m[pos].vId1 = 0;
      }
      else if(m[pos].vId2 == vID)
      {
        m[pos].vId2 = 0;
      }

      if(m[pos].vId1 == 0 && m[pos].vId2 == 0)
      {
        m.erase(pos);
      }
    }
    else
    {
      //assert(0);
    }
  }
  void detectAccident()
  {
    for(map<int, accAtPos>::iterator ii = m.begin(); ii != m.end(); ++ii)
    {
      if((*ii).second.vId1 != 0 && (*ii).second.vId2 != 0)
      {
#ifdef printAccidents
        //cout << "Accident happend at time=" << (*ii).second.time << " pos=" << (*ii).first << " for car1=" << (*ii).second.vId1 << " and car2=" << (*ii).second.vId2 << endl;
#endif
        //occuredAccidents++;
      }
    }
  }
  bool isAccidentInSeq(int seq)
  {
    for(map<int, accAtPos>::iterator ii = m.begin(); ii != m.end(); ++ii)
    {
      if((*ii).first / SEGMENT_LENGTH == seq)
      {
        if((*ii).second.vId1 != 0 && (*ii).second.vId2 != 0)
        {
            // cout << "acc found" << endl;
          return true;
        }
      }
    }
    return false;
  }
  size_t addFirst;
  size_t addSecond;
  size_t addNew;
  size_t addAlreadyFirst;
  size_t addAlreadySecond;
  size_t removeFirst;
  size_t removeSecond;
  size_t erase;
  size_t couldNotBe;
  size_t deleteCallCnt;
  size_t addCallCnt;

  size_t size() {
    return m.size();
  }
private:
  std::map<int, accAtPos> m;

};

struct LinearRoadTollRecord {

  int m_iType;
  int m_iVid;
  int m_iTime;
  int m_iEmit;
  int m_iSpd;
  int m_iToll;

  LinearRoadTollRecord() {
    m_iType = m_iVid = m_iTime = m_iEmit = m_iSpd = m_iToll = 0;
  }

  LinearRoadTollRecord(int type, int time, int vid, int emit, int spd, int toll) {
    m_iType = type;

    m_iTime = time;
    m_iVid = vid;
    m_iEmit = emit;
    m_iSpd = spd;
    m_iToll = toll;
  }
  LinearRoadTollRecord(const LinearRoadRecord &rec) {
    m_iType = rec.m_iType;
    m_iVid = rec.m_iVid;
    m_iTime = rec.m_iTime;
    /* m_iEmit = emit; */
    /* m_iSpd = spd; */
    /* m_iToll = toll; */
    m_iEmit = 0;
    m_iSpd = 0;
    m_iToll = 0;
  }

  std::string toString() {
    std::ostringstream oss;
    oss << "Type: " << m_iType << "\tTime: " << m_iTime << "\tEmit: " << m_iEmit
        << "\tSpd: " << m_iSpd << "\tToll: " << m_iToll;
    return oss.str();
  }
};

#endif // LRB_TOLL_RECORD_H
