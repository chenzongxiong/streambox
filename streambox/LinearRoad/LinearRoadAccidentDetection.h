#ifndef LRB_ACCIDENT_DETECTION_H
#define LRB_ACCIDENT_DETECTION_H

#include <map>


class AccidentDetection {
public:
  struct __attribute__((packed)) accAtPos
  {
    int vId1;
    int vId2;
    int time;
  };

  AccidentDetection()
  {
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

  void addStop(int vID, int pos, int time) {
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
          cout << "acc found" << endl;
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


private:
  std::map<int, accAtPos> m;  // pos
};


#endif // LRB_ACCIDENTDETECTION_H
