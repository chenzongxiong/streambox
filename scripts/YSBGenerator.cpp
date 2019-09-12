#include <iostream>
#include <cstring>
#include <omp.h>
#include <chrono>
#include <random>
#include <bitset>
#include <algorithm>
#include <atomic>
#include <array>
#include <assert.h>
#include <sstream>
#include <fstream>

using namespace std;

static const std::string events[] = {"view", "click", "purchase"};

struct __attribute__((packed)) record {
    uint8_t user_id[16];
    uint8_t page_id[16];
    uint8_t campaign_id[16];
    char ad_type[9];
    char event_type[9];
    int64_t current_ms;
    uint32_t ip;

    record(){
        event_type[0] = '-';//invalid record
        current_ms = 0;
        ip = 0;
    }

    record(const record& rhs)
    {
        memcpy(&user_id, &rhs.user_id, 16);
        memcpy(&page_id, &rhs.page_id, 16);
        memcpy(&campaign_id, &rhs.campaign_id, 16);
        memcpy(&event_type, &rhs.event_type, 9);
        memcpy(&ad_type, &rhs.ad_type, 9);
        current_ms = rhs.current_ms;
        ip = rhs.current_ms;
    }
};//size 78 bytes

void shuffle(record* array, size_t n)
{
    if (n > 1)
    {
        size_t i;
        for (i = 0; i < n - 1; i++)
        {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            record t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

void generate(record& data, size_t campaingOffset, uint64_t campaign_lsb, uint64_t campaign_msb, size_t event_id)
{
    event_id = event_id % 3;

    memcpy(data.campaign_id, &campaign_msb, 8);

    uint64_t campaign_lsbr = campaign_lsb + campaingOffset;
    memcpy(&data.campaign_id[8], &campaign_lsbr, 8);

    const char* str = events[event_id].c_str();
    //data.ad_type = "banner78";
    strcpy(&data.ad_type[0], "banner78");
    strcpy(&data.event_type[0], str);

    auto ts = std::chrono::system_clock::now().time_since_epoch();
    data.current_ms = std::chrono::duration_cast<std::chrono::milliseconds>(ts).count();

    data.ip = 0x01020304;
}


int main(int argc, char *argv[])
{
    cout << "usage filePath tupleCnt threadCnt campaingCnt" << endl;

    char* filePath = "file.bin";
    size_t tupleCnt = 1000;
    size_t threadCnt = 1;
    size_t campaingCnt = 1000;

    if(argc == 5)
    {
        filePath = argv[1];
        tupleCnt = atoi(argv[2]);
        threadCnt = atoi(argv[3]);
        campaingCnt = atoi(argv[4]);
        cout << "using paramter filepath=" << filePath << " tupleCnt=" << tupleCnt << " threadCnt=" << threadCnt << " campaingCnt=" << campaingCnt << endl;
    }
    else
    {
        cout << "not enough parameter" << endl;
    }

    std::vector<string> fileNames;
    for(size_t i = 0; i < threadCnt; i++)
    {
        std::string path(filePath);
        path = path.substr(0,path.find("."));
        path.append(to_string(i));
        path.append(".bin");
        fileNames.push_back(path);
        cout << "fileName " << i << "=" << path << endl;
    }

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<size_t> diss(0, SIZE_MAX);

    size_t randomCnt = tupleCnt/10;
    size_t* randomNumbers = new size_t[randomCnt];
    std::uniform_int_distribution<size_t> disi(0, campaingCnt);
    for(size_t i = 0; i < randomCnt; i++)
        randomNumbers[i] = disi(gen);

    record** recs = new record*[threadCnt];
    for(size_t i = 0; i < threadCnt; i++)
    {
        recs[i] = new record[tupleCnt];
    }


    uint64_t campaign_lsb, campaign_msb;
    auto uuid = diss(gen);
    uint8_t* uuid_ptr = reinterpret_cast<uint8_t*>(&uuid);
    memcpy(&campaign_msb, uuid_ptr, 8);
    memcpy(&campaign_lsb, uuid_ptr + 8, 8);
    campaign_lsb &= 0xffffffff00000000;

#pragma omp parallel num_threads(threadCnt)
    {
#pragma omp for
        for(size_t th = 0; th < threadCnt; th++)
        {
            for(size_t i = 0; i < tupleCnt; i++)
            {
                generate(recs[omp_get_thread_num()][i], /**campaingOffset*/ randomNumbers[i%randomCnt], campaign_lsb, campaign_msb, /**eventID*/ i);
                //cout << "thread=" << omp_get_thread_num() << " put in " << i << " value=" << recs[omp_get_thread_num()][i].event_type << endl;
            }

            shuffle(recs[omp_get_thread_num()], tupleCnt);

            //write to file
            cout << "writing file" << fileNames[omp_get_thread_num()] << endl;
            std::ofstream ofp(fileNames[omp_get_thread_num()], std::ios::out | std::ios::binary);
            ofp.write(reinterpret_cast<const char*>(recs[omp_get_thread_num()]), tupleCnt * sizeof(record));
            ofp.close();

            //DBG read first 10
            //    record* inRecs = new record[tupleCnt];
            //    std::ifstream ifp(fileNames[omp_get_thread_num()], std::ios::in | std::ios::binary);
            //    ifp.read(reinterpret_cast<char*>(inRecs), tupleCnt * sizeof(record));
            //    for(size_t i = 0; i < tupleCnt; i++)
            //    {
            //      cout << "eventtype " << i << "=" << inRecs[i].event_type << endl;
            //    }
            //    ifp.close();
        }

    }
}
